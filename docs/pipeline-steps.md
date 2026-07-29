# Pipeline steps

The BEguided pipeline is a sequence of six steps. Steps 1–3 always run. Steps 4–5 only run if `BEguidedOFF` is installed and `search_offtargets: true` is set. Step 6 is the final merge.

```
inputs (TSV + FASTA + GTF + base editor profiles)
        │
        ▼
┌─────────────────────────────────────┐
│ Step 1: guidefinder                 │  ← scan for PAM sites, build guides
│   in : genome FASTA, targets TSV,   │
│        base editor profiles TSV     │
│   out: step1_guidefinder.tsv        │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 2: bystander                   │  ← predict bystander edits in CDS
│   in : step1 + GTF annotation       │
│   out: step2_bystander.tsv          │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 3: filter                      │  ← CRISPRater, GC, homopolymer,
│   in : step2                        │     trinuc context score
│   out: step3_filtered.tsv           │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 4: off-target search           │  ← BEguidedOFF + CCTOP scoring
│   in : step3                        │     (optional, GPU-accelerated)
│   out: step4_offtarget.tsv          │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 5: off-target filter           │  ← drop guides with too many / too
│   in : step4                        │     concerning off-targets
│   out: (merged into step 6)         │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 6: merge                       │  ← combine step3 + step5 into
│   in : step3, step4                 │     the final TSV
│   out: beguided_results.tsv         │
└─────────────────────────────────────┘
```

## Step 1 — Guide generation (`step1_guidefinder.tsv`)

For every target row in the input TSV and every base editor profile, this step:

1. Loads the genome chunk around the target position.
2. Searches both strands for valid PAM sites matching any profile's PAM (IUPAC-aware).
3. For each PAM hit, takes the 20-nt (or profile-specific) protospacer on the PAM-proximal side.
4. Inserts the desired base at the editing position implied by the profile's editing window, producing a hypothetical edited sequence.
5. Filters to guides where the desired edit lies inside the profile's *high-efficiency* window (else the guide is generated but the editing efficiency would be too low).
6. Computes the CRISPRater on-target activity score (`guide_score`) for 20-nt guides.
7. Computes `gc_content`, `homopolymer`, `trinuc_edit_score` for every guide length.

Output row count: targets × guides per target × sense+antisense variants of each editor × matching PAMs.

## Step 2 — Bystander analysis (`step2_bystander.tsv`)

For every guide produced in step 1, this step:

1. Determines the protospacer sequence and the editing window implied by the base editor profile.
2. Slides every other editable base within the window (per the profile's editing chemistry) and checks the codon context using the supplied GTF annotation.
3. Classifies each bystander as silent (same amino acid), missense (different amino acid), or nonsense (premature stop).
4. Attaches per-bystander metadata: position (1-based from the 5' end of the guide), original codon, edited codon, original amino acid, edited amino acid, classification, transcript / exon context.

Step 2 only runs when the genome FASTA and GTF annotation are both available. The `keep_bystanders` config knob controls which guides step 3 keeps:
- `all` — keep every guide regardless of bystander status
- `no_bystanders` — keep only guides with zero bystanders
- `silent` — keep guides whose bystanders are all silent (or absent)

## Step 3 — Filtering (`step3_filtered.tsv`)

Applies the filters configured under the pipeline config:

1. `guide_score_threshold` — drop guides below the threshold
2. `gc_content_filter` — drop by GC bucket
3. `filter_homopolymers` — drop guides with ≥4-mer homopolymer runs
4. `trinuc_edit_score_min` — drop ABE guides below the trinucleotide context threshold (CBE guides unaffected)

Every filter that drops guides is logged at INFO so you can see why a guide disappeared.

## Step 4 — Off-target search (`step4_offtarget.tsv`)

**Optional.** Skipped if `BEguidedOFF` is not installed or `search_offtargets: false`.

For every guide that survives step 3, this step:

1. Encodes the guide against the pre-built BEguidedOFF `.npz` genome.
2. Hamming-distance search up to `mismatch_tolerance` mismatches in the 20-nt guide region.
3. Scores every off-target with CCTOP (penalty model that weighs mismatches by PAM proximity).
4. Keeps the top-1 (worst) off-target per guide to minimise output size.

## Step 5 — Off-target filter

For every guide with at least one off-target in step 4, applies:

- `remove_guides_score` — drop guides whose worst off-target score is below this value (i.e. too concerning)
- `remove_guides_n_offtargets` — drop guides whose count of concerning off-targets exceeds this number (if set)

Step 5 doesn't write a separate file; it runs inline and the resulting keep/discard decisions feed into step 6.

## Step 6 — Merge (`beguided_results.tsv`)

Joins step 3 (filtered guides) with step 5 (off-target decisions) and writes the final TSV. Columns:

- Everything from step 3 (guide metadata, bystander info, sequence properties)
- The five off-target summary columns (`worst_offtarget_score`, `offtarget_count`, `worst_offtarget_chr`, `worst_offtarget_pos`, `worst_offtarget_sequence`) when applicable; `NaN` if off-target search was skipped
- All extra columns from the input TSV (`passthrough_extras`)
