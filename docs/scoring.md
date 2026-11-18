# Scoring semantics

BEguided produces four scores per guide. This page explains what each one means, how it's computed, and how to use it.

## `guide_score` — CRISPRater on-target activity

- **Range:** 0–1, higher = more active guide.
- **When:** step 1.
- **How:** Re-implementation of the CRISPRater gradient-boosted model (originally from CCTOP standalone). Predicts the likelihood that the guide drives strong editing in a typical target context.
- **Constraint:** Only computed for **20-nt guides**. Every bundled base editor profile in `configs/base_editor_profiles.tsv` uses 20-nt guides, so the column is always populated. If you add a custom profile with a different length, the column is left empty with a DEBUG log line.
- **Filter:** `guide_score_threshold` (default `0.5`) drops guides below this value in step 3. Set to `0` to disable.

## `gc_content` — protospacer GC

- **Range:** 0–100 (percent).
- **When:** step 1, every guide length.
- **How:** `(G + C) / length * 100` of the 20-nt (or profile-specific) protospacer.
- **Filter:** `gc_content_filter` (`off` / `high` / `medium`). `high` keeps 40–60%, `medium` keeps 30–70%. The deprecated names `good` (= `high`) and `acceptable` (= `medium`) are still accepted. Guides with a missing GC value (shouldn't happen with DNA input) are always kept.

## `homopolymer` — ≥4-mer run flag

- **Type:** boolean.
- **When:** step 1, every guide length.
- **How:** `true` if the guide contains any run of 4 or more identical bases (`AAAA`, `TTTT`, `CCCC`, `GGGG`). The flag covers both Pol III termination signals (A/T runs) and potential primer-dimer artefacts (C/G runs).
- **Filter:** `filter_homopolymers` (default `false`). Set to `true` to drop flagged guides in step 3.

## `trinuc_edit_score` — Marquart et al. 2021 ABE8e context

- **Range:** 0–1, higher = more likely to edit. `CAC` = 1.0, `AAA` ≈ 0.02.
- **When:** step 1, every ABE guide. CBE guides are always `NaN`.
- **How:** Lookup in the bundled, normalised Marquart et al. 2021 (ABE8e) `X-A-Y` trinucleotide table. The lookup uses the mRNA-sense context: for `-` strand guides the genomic 3-mer is reverse-complemented before lookup.
- **Caveat:** The bundled table is derived from ABE8e spCas9 data and is applied to every ABE profile as an approximation. Actual editing efficiencies will differ between ABE8e, ABE9, ABE-Umax, etc. CBE profiles always produce `NaN` because the table is ABE-only. There is no config override — the bundled table is always used.
- **Filter:** `trinuc_edit_score_min` (default `null`) drops ABE guides below this threshold. CBE guides unaffected.

## `worst_offtarget_score` — CCTOP penalty score

- **Range:** 0 (perfect match) to large positive (unconcerning). No off-targets = blank (`NaN`).
- **When:** step 4, only when `BEguidedOFF` is installed and `search_offtargets: true`.
- **How:** CCTOP penalty model. Mismatches near the PAM are penalised more than mismatches at the PAM-distal end of the guide. The score is unbiased across PAM families.
- **Filter:** `remove_guides_score` (default `5.0`) drops guides whose worst off-target score is below this value (i.e. too concerning). The semantics are "lower = worse", so a threshold of `5.0` keeps guides whose worst off-target scores `5.0` or higher. Guides with no off-targets (blank score) always pass.

## `offtarget_count` — concerning off-target count

- **Range:** non-negative integer. 0 = no concerning off-targets.
- **When:** step 4.
- **How:** Number of off-targets with CCTOP score ≤ `count_threshold_score` (default `15.0`).
- **Filter:** `remove_guides_n_offtargets` (default `100`; `null` = no limit) drops guides whose count exceeds this number.

## Score ranges cheat sheet

| Score | Better value | Typical cutoff |
|---|---|---|
| `guide_score` | higher | `>= 0.5` |
| `gc_content` | around 50% | `40–60%` or `30–70%` |
| `trinuc_edit_score` | higher | `>= 0.5` (ABE only) |
| `worst_offtarget_score` | higher | `>= 5.0` |
| `offtarget_count` | lower | `<= 100` (if set) |

> **Putting it together:** for an ABE guide against a typical target, a high-quality result is `guide_score >= 0.75`, `gc_content in [40, 60]`, `homopolymer == false`, `trinuc_edit_score >= 0.5`, and `worst_offtarget_score >= 10`. Use the dynamic filters in the web UI to explore the trade-offs without re-running the pipeline.
