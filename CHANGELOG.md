# Changelog

All notable changes to verkko-fillet will be documented in this file.

The project follows [Semantic Versioning](https://semver.org/).

## [0.1.28] - 2026-10-09

### Fixed
- telomere_analysis.sh were fixed to locate the correct directory
- remove verkko_hap

## [0.1.27] - 2026-10-01

### Added
- new function gapCleaning, which wrapper of map_rDNA, find_gap_in_rDNA, rDNA_gap_cleaning
- new function gapCleaning, which wrapper of map_rDNA, find_gap_in_rDNA, rDNA_gap_cleaning

### Fixed
- readChr function accept haplotype1 and haplotype2 with str, and check if contig names contain the haplotype name, rather parsing the contig names.
- n50Plot function uses hap column, instead of hap_verkko column from obj.stats
- functions and scripts that are associated with internal telomere detection are fixed. The paths of the execs were wrong
- functions and scripts that are associated with internal telomere detection are fixed. The paths of the execs were wrong
- find_telomere.sh can find the correct vgp scripts
- read_Verkko can ignore assembly.paths.tsv
- make_verkko_fillet_dir.sh bug fix. make verkko-fillet directory after remove existing one

## [0.1.26] - 2026-09-28

### Fixed
- bugs
- run_shell function will show exact error with raise
- removeRDNA.sh fix typo when checking mash
- readNode and read_Verkko function will ignore absence of the assembly.colors.csv for non-HiC or non-trio verkko assemblies
- chromosome assignment script accept new parameter min_Length

## 0.1.26

### Added

- Initial changelog tracking for verkko-fillet.
