# Alternative mutation sources

Dataset preparation — downloading public mutation databases and converting
them into BEguided-ready target lists (`chromosome`, `position`,
`target_base`, `desired_base`) — lives in the companion repository
**BEguided-datasets**, a sibling checkout of this repo.

Supported sources: COSMIC (BEguided's reference inputs), OncoDB, Cancer
Hotspots, DORM, ClinVar, and cBioPortal.

Everything related lives there:

| | |
|---|---|
| Converter scripts (`build_<source>_targets.py`, `protein_to_genomic.py`) | `../BEguided-datasets/scripts/` |
| Per-source documentation | `../BEguided-datasets/docs/alternative-mutation-sources.md` |
| Full-run configs (per source × ABE/CBE) | `../BEguided-datasets/configs/` |
| Smoke configs (first-5-rows sanity runs) | `../BEguided-datasets/smoke/` |
| Dataset version manifest | `../BEguided-datasets/datasets.yaml` |

BEguided itself only consumes the generated TSVs through the standard input
contract — see [input-output.md](input-output.md). The pipeline has no
dependency on the dataset repo; point any config's `input.targets` at a
generated list and run.
