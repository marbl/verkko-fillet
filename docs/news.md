(News)=

# 📣 News

<!-- marker: after prelude -->

## 🎉 v0.1.23 released

The latest release of `verkko-fillet` is out, with a new automated preprocessing CLI, gzipped FASTA support across the chromosome-renaming pipeline, mashmap caching, plotting improvements (PDF by default, `force` overwrite, tighter layouts), and several QoL fixes.

- 📄 Full changelog: [v0.1.23](release-notes/0.1.23.md)
- 🆕 New tutorial: [Automatic Preprocessing](tutorials/basics/auto_preprocessing.md) — run the standard QC + chromosome-assignment + telomere pipeline end-to-end with a single command.

We strongly recommend upgrading to the latest version for refining and cleaning Verkko assemblies.

---

<!-- marker: before old news -->

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