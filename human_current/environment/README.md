# Environment notes

This code collection spans several production lanes and does not represent one
uniformly validated software environment.

## Legacy v5r2 environment

`legacy_v5r2/environment.yml` and `legacy_v5r2/requirements.txt` are the
retained environment declarations from the earlier public release. They specify
Python 3.11.9, R 4.3.3 and the versions listed in those files.

## Formal21 regional-model environment

The formal regional-model snapshot was recorded on Ubuntu 22.04.5 LTS with
Python 3.10.12 and R 4.5.0. `formal21/formal21_versions.tsv` contains the
relevant package versions selected from that production snapshot. The original
machine-specific library paths are not part of this release.

## Other current lanes

The fixed-H external scorer requires NumPy, SciPy and PyTorch with the approved
CUDA route. Current figure producers additionally use packages named in the
main README, including PyMuPDF, pypdf, Pillow, Matplotlib, PyArrow and OpenPyXL.
R figure and local-contrast producers use packages including lme4, pbkrtest,
ComplexHeatmap and circlize. PDF preview generation uses Poppler, and several
figures use installed DejaVu Sans or Liberation Sans fonts.

Exact package snapshots were not retained for every recent figure-composition
lane. The legacy and formal21 files must therefore be read as lane-specific
records, not as a single lock file covering every script.
