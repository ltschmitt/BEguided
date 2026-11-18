# Input and output formats

## Input: target file

The target input file is a tab-separated file (`.tsv`) with four required columns:

| Column | Description |
|---|---|
| `chromosome` | Chromosome name (e.g. `chr1`, `12`, `X`) |
| `position` | 1-based genomic coordinate of the target mutation |
| `target_base` | Base currently present at the target (`A`, `T`, `C`, `G`) |
| `desired_base` | Desired base after editing (`A`, `T`, `C`, `G`) |

You may add any number of extra columns (`gene`, `mutation_occurrence`, `genomic_mutation_id`, etc.); they are carried through every pipeline step into the final output (`passthrough_extras: true` by default).

Example:

```tsv
chromosome	position	target_base	desired_base
chr1	50	T	C
chr1	75	A	G
chr2	30	T	C
chr2	80	A	G
```

The web UI accepts `.tsv`, `.csv`, and `.xlsx/.xls`. The CLI only reads `.tsv`. CSV / Excel inputs are normalised to TSV on the server side.

### Sequence input mode (web UI only)

The web UI also accepts pasted sequences with a bracket marker for the desired edit, e.g. `ATAG[T>C]GACGTACGT…`. The pasted sequence serves as the reference, so guide scoring runs on the user's sequence; bystander analysis, which requires a genomic transcript annotation, is not performed and results carry no genomic coordinates. See [`web-deployment.md`](web-deployment.md) for details.

### Single-target input (web UI only)

In the web UI's *single* mode, one target is entered manually as chromosome, 1-based position, target base and desired base. It is processed exactly like a one-row spreadsheet upload (`POST /single/parse` normalises it to the same payload).

## Output: file layout

When you run the full pipeline with `output_prefix: "./results"`, files are written directly into that directory:

```
results/
├── step1_guidefinder.tsv     ← step 1: candidate guide RNAs
├── step2_bystander.tsv       ← step 2: bystander analysis
├── step3_filtered.tsv        ← step 3: filtered guides
├── step4_offtarget_input.tsv ← step 4: BEguidedOFF input (optional)
├── step4_offtarget.tsv       ← step 4: BEguidedOFF off-targets (optional)
└── beguided_results.tsv      ← step 5: final merged output
```

The `keep_intermediates` config key is accepted, but per-step cleanup is not implemented yet — the intermediate step files are always kept on disk.

## Output: column reference (final `beguided_results.tsv`)

The final TSV contains the merged guide + off-target summary. The exact column set depends on your config; common columns:

### Core (always present)

| Column | Description |
|---|---|
| `chromosome` | From input |
| `position` | From input |
| `target_base` | From input |
| `desired_base` | From input |
| `base_editor` | The editor that produced this guide (sense + antisense variants of each profile) |
| `guide_rna_sequence` | The 20-nt protospacer (always 20 bp for the bundled profiles) |
| `pam_sequence` | The PAM detected adjacent to the guide |
| `editing_position` | 1-based position of the edited base within the guide (from 5') |
| `strand_used` | `+` or `-` (which strand carries the PAM) |

### Sequence-property columns (always present)

| Column | Description |
|---|---|
| `guide_score` | CRISPRater on-target activity score (0–1, higher = better). **Only computed for 20-nt guides**; non-20-nt profiles are skipped with a warning and the column is empty. |
| `gc_content` | GC percentage of the whole guide (0–100) |
| `homopolymer` | `true` if the guide contains any run of ≥4 identical bases |
| `trinuc_edit_score` | Marquart et al. 2021 normalised ABE8e proportion for the local `X-A-Y` context (0–1). `NaN` for CBE guides (table is ABE-only). The `-` strand is reverse-complemented before lookup so the value reflects the mRNA-sense context the paper measured. |

### Off-target columns (only when `search_offtargets: true` and `BEguidedOFF` is installed)

| Column | Description |
|---|---|
| `worst_offtarget_score` | Lowest CCTOP score (most dangerous; blank/`NaN` if no off-targets) |
| `offtarget_count` | Number of off-targets with CCTOP score ≤ `count_threshold_score` |
| `worst_offtarget_chr` | Chromosome of the worst off-target |
| `worst_offtarget_pos` | 1-based position of the worst off-target |
| `worst_offtarget_sequence` | Sequence of the worst off-target |

CCTOP score semantics: **lower = more concerning** (mismatches near PAM). Perfect match = 0.0. No off-targets = blank (`NaN`; always passes `worst_offtarget_score >= X` filters).

### Bystander columns (from step 2)

Step 2 attaches one row per *guide × base_editor* (since each editor has its own bystander window). Columns include the predicted bystander positions, the resulting codon changes, the amino-acid substitutions and the `bystander_effect` classification: `none`, `synonymous`, or `nonsynonymous` (any amino-acid change, including nonsense/stop-gain), plus the labels `multi-transcript`, `editing_position_non_coding`, `no_cds_overlap`, `bystander_outside_cds` and `cds_reconstruction_failed` for cases where no codon-level consequence can be assigned. See [`pipeline-steps.md`](pipeline-steps.md) for the full label reference.

### Passthrough columns

Any extra column from your input TSV is preserved through every step into the final file.
