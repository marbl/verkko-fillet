(News)=

# 📣 News

<!-- marker: after prelude -->

## 🎉 v0.1.26 released

This release includes several bug fixes and improvements based on recently reported issues. Thank you to everyone who reported issues and provided feedback!

Bug Fixes

- **Improved Bash script execution**
  - Fixed an issue where Bash scripts were executed using `sh` from Python subprocess calls.
  - All Bash-based scripts are now explicitly executed with `bash`.
  - This should resolve script execution and symbolic-link creation errors on systems where `/bin/sh` points to a non-Bash shell, such as `dash`.
  - Related to #11 and #13.

- **Support for assemblies without `assembly.colors.tsv`**
  - verkko-fillet can now handle assemblies where `assembly.colors.tsv` is not available.
  - This allows unphased Verkko assemblies to be loaded without raising an error.
  - Fixes #14.

Improvements

- **Configurable minimum contig length in `chrAssign`**
  - Added the `min_Length` parameter to `chrAssign`.
  - The default remains `5,000,000` bp to preserve the previous behavior.
  - Users can now lower the threshold for assemblies containing smaller chromosomes, such as bird dot chromosomes.
  - Addresses #12.

We recommend updating to the latest version and rerunning any steps affected by these issues.

Thank you again for your feedback and contributions to improving verkko-fillet!

## 🎉 v0.1.25 released

This release adds new command-line utilities for contaminant screening and rDNA gap cleaning, plus clearer pre-run documentation for reproducible Verkko setup.

- 🧪 New scripts: `screen-assembly.sh` and `rDNA_gap_cleaning.sh` for standardized assembly screening and rDNA-bounded gap replacement.
- 📚 Documentation updates: improved scientific guidance in the new pre-run tutorial pages for input checks and Verkko execution modes.
- ⚙️ Reliability improvements: stronger input validation, processing safeguards, and progress reporting in the new utilities.

📄 Full changelog: [v0.1.25](release-notes/0.1.25.md)

🆕 Pre-run tutorials:
- [Check Input](tutorials/before_start/check_input.md)
- [Run Verkko](tutorials/before_start/run_verkko.md)

<!-- marker: before old news -->

## 🎉 v0.1.24 released

A small maintenance release with debugging and reliability improvements to `getChrNames.sh`:
- The `neighborhood` executable bundled alongside the script is now picked up automatically — no manual install required.

📄 Full changelog: [v0.1.24](release-notes/0.1.24.md)

## 🎉 v0.1.23 released

The latest release of `verkko-fillet` is out, with a new automated preprocessing CLI, gzipped FASTA support across the chromosome-renaming pipeline, mashmap caching, plotting improvements (PDF by default, `force` overwrite, tighter layouts), and several QoL fixes.

- 📄 Full changelog: [v0.1.23](release-notes/0.1.23.md)
- 🆕 New tutorial: [Automatic Preprocessing](tutorials/basics/auto_preprocessing.md) — run the standard QC + chromosome-assignment + telomere pipeline end-to-end with a single command.

We strongly recommend upgrading to the latest version for refining and cleaning Verkko assemblies.



---

## 📚 Previous highlights

### v0.1.21
This release shipped numerous updates, new features, and bug fixes. See the [v0.1.21 release notes](release-notes/0.1.21.md).

### Tutorial: Recovering T2T contigs
A comprehensive [tutorial on recovering Telomere-to-Telomere (T2T) assemblies](tutorials/basics/telo.ipynb) is available. It covers detecting internal telomeres, reconnecting broken contigs using chromosome-assignment information, and linking small nodes to main contigs via graph alignment.

### v0.1.19
Bug fixes and refinements. See the [v0.1.19 release notes](release-notes/0.1.19.md).

### v0.1.18 — major update
Numerous new functions and `FilletObj` attributes. See the [v0.1.18 release notes](release-notes/0.1.18.md).

> **Note:** v0.1.18 introduced new Python dependencies, including `scikit-learn`.

### `verkko-fillet` on PyPI — *2025-01-21*
`verkko-fillet` is now installable via `pip`:

```bash
pip install verkkofillet
```