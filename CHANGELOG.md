# CHANGELOG

## [unreleased]
- Handle SGL (single end BNDs) - ignore unless large insertions (#70)
- Lint and format code with ruff, add an action (#76)
- Keep a changelog action (#76)
- Run smoke test action on pull requests (#76)
- Trim `chr` prefix from blacklist file (#77)
- Provide an alias exclusionlist for blacklist, being deprecated (#77)

## [0.9.3]
- Update release automation build environment (#67)

## [0.9.2]
- Remove any chr prefix from input chromosome contig name (#65)

## [0.9]
- Blacklist (excluded regions list) also hides probes

## [0.8]
- Changed coverage probe height calculation to a log2-ratio
- Aberration `gain` is no longer set to `true` for deletions (giving "loss" in vcf2cytosure)
- Changed coverage average calculation to also exclude 0-coverage bins. Note that N-sequence bins already get no probes.
- Add Dockerfile
