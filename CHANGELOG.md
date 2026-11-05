# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [2.0.0] - 2026-11-05

### Breaking
- RNA expression/classification/QC moved from the sample document to the `rna_expression`, `rna_classification` and `rna_qc` collections.
- Module-level `load_*` / `update_case` functions replaced by `IngestionService`.
- Translocation import now keeps fusion-relevant DEL/DUP SVs (previously dropped).

### Added
- Class-based ingestion: `IngestionService`, `DnaParser`, `RnaParser`, `DependentWriter`, `Repos`.
- Rollback of partial imports on error (`cleanup_on_error`).
- `ingest_status` and `data_counts` on sample documents.
- `--dev` flag for the `yaml` command; active database is logged.
- SNV filter: skip `FAIL_NVAF`, `FAIL_LONGDEL`, `FAIL_PON_*`.
- Translocations: BND `feature_fusion` mate merging, `transcript_ablation` (BND only), PASS-only, MANE annotation.

### Fixed
- 1000 Genomes `GMAF` parsed per allele.
- Case/control `null` handling.

## [1.0.0] - 2026-10-02

Legacy Python CLI replacing the four Perl import scripts.

- `load` (flags) and `yaml` subcommands; DNA and RNA samples.
- DNA: SNVs/indels, CNVs, biomarkers, low coverage, coverage, translocations. RNA: fusions, expression, classification, QC.
- `--increment` case ID suffixing and `--update` for existing cases.
- CSQ reduction to selected transcript; population AF summary (gnomAD, ExAC, 1000G).
- Database `coyote3`; `cnvs_wgs` renamed to `cnvs`.
- `migrate_variants_idref.py` and `DATA_INPUT_SPEC.md`.

[2.0.0]: https://github.com/SMD-Bioinformatics-Lund/coyote_cli/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/SMD-Bioinformatics-Lund/coyote_cli/releases/tag/v1.0.0
