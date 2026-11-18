# Pipeline steps

The BEguided pipeline is a sequence of five steps. Steps 1–3 always run. Steps 4–5 only run if `search_offtargets: true` is set, which requires `BEguidedOFF` to be installed (the run aborts with an error otherwise). Step 5 performs the off-target filtering and the final merge inline.

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
│   in : step3                        │     (optional; CPU "fast" mode or
│   out: step4_offtarget.tsv          │     GPU "extensive" mode)
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│ Step 5: off-target filter & merge   │  ← drop guides with too many / too
│   in : step3, step4                 │     concerning off-targets, then
│   out: beguided_results.tsv         │     combine with step3 into final TSV
└─────────────────────────────────────┘
```

## Step 1 — Guide generation (`step1_guidefinder.tsv`)

For every target row in the input TSV and every base editor profile, this step:

1. Loads the genome chunk around the target position.
2. Searches both strands for valid PAM sites matching any profile's PAM (IUPAC-aware).
3. For each PAM hit, takes the 20-nt (or profile-specific) protospacer on the PAM-proximal side.
4. Inserts the desired base at the editing position implied by the profile's editing window, producing a hypothetical edited sequence.
5. Keeps only guides whose desired edit lies inside the profile's *high-activity* window (`high_edit_start`–`high_edit_end`); a target outside this window is not reported as editable by that profile.
6. Computes the CRISPRater on-target activity score (`guide_score`) for 20-nt guides.
7. Computes `gc_content`, `homopolymer`, `trinuc_edit_score` for every guide length.

Output row count: targets × guides per target × sense+antisense variants of each editor × matching PAMs.

## Step 2 — Bystander analysis (`step2_bystander.tsv`)

For every guide produced in step 1, this step:

1. Determines the protospacer sequence and the *general* editing window (`general_edit_start`–`general_edit_end`) implied by the base editor profile.
2. Slides every other editable base within the window (per the profile's editing chemistry) and checks the codon context using the supplied GTF annotation.
3. Classifies each bystander by its protein-level consequence: synonymous (same amino acid) or nonsynonymous (any amino-acid change). Missense and nonsense (stop-gain) are collapsed into this single class — the finer distinction is not maintained.
4. Attaches per-bystander metadata: position (1-based from the 5' end of the guide), original codon, edited codon, original amino acid, edited amino acid, classification, transcript / exon context.

The per-guide `bystander_effect` column carries one of the following labels:

| Label | Meaning |
|---|---|
| `none` | No editable bystander positions |
| `synonymous` | All bystander edits are synonymous |
| `nonsynonymous` | At least one bystander changes the amino acid (missense and nonsense/stop-gain collapsed) |
| `multi-transcript` | The guide overlaps CDS annotations of several transcripts, so a single per-transcript consequence cannot be assigned |
| `editing_position_non_coding` | The target base lies outside an annotated CDS |
| `no_cds_overlap` | The guide overlaps no annotated CDS |
| `bystander_outside_cds` | A bystander falls outside the CDS exon the guide overlaps |
| `cds_reconstruction_failed` | A codon spans an exon junction and could not be reconstructed from the available data |

`none`, `synonymous` and `nonsynonymous` are the coding consequences; the remaining labels mark cases where no codon-level consequence could be assigned (non-coding or ambiguous context).

Step 2 only runs when the genome FASTA and GTF annotation are both available. The `keep_bystanders` config knob controls which guides step 3 keeps:
- `all` — keep every guide regardless of bystander status
- `none` — keep only guides with zero bystanders
- `synonymous` — keep guides whose bystanders are all synonymous (or absent)

## Step 3 — Filtering (`step3_filtered.tsv`)

Applies the filters configured under the pipeline config:

1. `guide_score_threshold` — drop guides below the threshold
2. `gc_content_filter` — drop by GC bucket (`high` = 40–60%, `medium` = 30–70%)
3. `filter_homopolymers` — drop guides with ≥4-mer homopolymer runs
4. `trinuc_edit_score_min` — drop ABE guides below the trinucleotide context threshold (CBE guides unaffected)

Every filter that drops guides is logged at INFO so you can see why a guide disappeared.

## Step 4 — Off-target search (`step4_offtarget.tsv`)

**Optional.** Runs only when `search_offtargets: true` — and fails with a clear error if `BEguidedOFF` is not installed. Skipped entirely when `search_offtargets: false`.

For every guide that survives step 3, this step:

1. Prepares the BEguidedOFF input with on-target exclusion positions (`step4_offtarget_input.tsv`).
2. Searches the reference genome with BEguidedOFF:
   - `extensive` mode (GPU): Hamming-distance search up to `offtargets.max_mismatches` mismatches in the 20-nt guide region against the pre-encoded `.npz` genome.
   - `fast` mode (CPU): FM-index search that enumerates exactly the sites with CCTOP score ≤ the search cutoff (derived from `remove_guides_score` / `count_threshold_score`; no GPU needed).
3. Scores every off-target with CCTOP (penalty model that weighs mismatches by PAM proximity).
4. Keeps the top-1 (worst) off-target per guide to minimise output size.

## Step 5 — Off-target filter & merge (`beguided_results.tsv`)

For every guide with at least one off-target in step 4, applies:

- `remove_guides_score` — drop guides whose worst off-target score is below this value (i.e. too concerning)
- `remove_guides_n_offtargets` — drop guides whose count of concerning off-targets exceeds this number (if set)

The filtering and the merge run inline in this step; no separate file is written for the filter. The step joins step 3 (filtered guides) with the off-target decisions and writes the final TSV. Columns:

- Everything from step 3 (guide metadata, bystander info, sequence properties)
- The five off-target summary columns (`worst_offtarget_score`, `offtarget_count`, `worst_offtarget_chr`, `worst_offtarget_pos`, `worst_offtarget_sequence`) when applicable; `NaN` if off-target search was skipped
- All extra columns from the input TSV (`passthrough_extras`)
