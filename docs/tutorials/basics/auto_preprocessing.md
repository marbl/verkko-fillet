# Automatic Basic Preprocessing

`automatic_script.py` is a command-line driver that runs the standard
verkko-fillet preprocessing and QC pipeline end-to-end on a finished Verkko
assembly. It is intended for users who do not want to step through the
notebook tutorials interactively and just want the figures, stats, and
chromosome-assignment files generated in one go.

The script is shipped with the package at:

```
verkko-fillet/src/verkkofillet/bin/automatic_script.py
```

## What it does

Given a Verkko output directory, a reference FASTA, and a chromosome name
mapping file, the script:

1. Loads the Verkko assembly into a `VerkkoFillet` object (`vf.pp.read_Verkko`)
2. Detects T2T contigs (`vf.tl.getT2T`)
3. Assigns reference chromosomes via mashmap (`vf.tl.chrAssign`)
4. Reads the chromosome mapping (`vf.pp.readChr`)
5. Detects broken contigs (`vf.pp.detectBrokenContigs`)
6. Generates QC plots:
   - `vf.pl.showMashmapOri` — chromosome coverage from mashmap
   - `vf.pl.completePlot` — completeness per chromosome/haplotype
   - `vf.pl.contigLenPlot` — contig length per haplotype
   - `vf.pl.contigPlot` — T2T status heatmap
   - `vf.pl.n50Plot` — N50 line plot
7. Detects internal telomeres and computes telomere percentages
   (`vf.tl.detect_internal_telomere`, `vf.pp.find_intra_telo`,
   `vf.pl.percTel`)

All outputs are written under the Verkko-fillet working directory
(`obj.verkko_fillet_dir`):

| Output                         | Location                                  |
| ------------------------------ | ----------------------------------------- |
| Figures (PDF)                  | `<verkko_fillet_dir>/figs/`               |
| Stats tables                   | `<verkko_fillet_dir>/stats/`              |
| Chromosome assignment files    | `<verkko_fillet_dir>/chromosome_assignment/` |

## Requirements

* `verkkofillet` installed (and importable as `verkkofillet`)
* External tools used by the underlying functions: `mashmap`, `samtools`,
  `bgzip`, `seqtk` (as required by `chrAssign` / FASTA renaming)
* A completed Verkko assembly directory

## Arguments

| Argument                    | Required | Default | Description                                                                 |
| --------------------------- | -------- | ------- | --------------------------------------------------------------------------- |
| `--verkkoDir`               | yes      | —       | Path to the Verkko output directory.                                        |
| `--ref`                     | yes      | —       | Reference FASTA file used for chromosome assignment.                        |
| `--map_file`                | yes      | —       | Tab-separated mapping from contig names in the reference FASTA to the desired chromosome names. |
| `--internal_tel_threshold`  | no       | `15000` | Distance from contig end (bp) used to flag internal telomeres.              |

### `--map_file` format

A two-column, tab-separated file:

```
<contig_name_in_fasta>\t<desired_chromosome_name>
```

Examples:

```
# contig is already named as desired
chr1        chr1
chr2        chr2     


# rename contig GCF00001 to chr1
GCF00001    chr1        
GCF00002    chr2
```

## Usage
```bash
python /path/to/verkko-fillet/src/verkkofillet/bin/automatic_script.py \
    --verkkoDir /path/to/verkko_output \
    --ref       /path/to/reference.fasta \
    --map_file  /path/to/chrMap.tsv \
    --internal_tel_threshold 15000
```

## Example

```bash
python automatic_script.py \
    --verkkoDir ./verkko_asm \
    --ref       ./ref/CHM13v2.fa \
    --map_file  ./ref/chrMap.tsv
```

After the run completes, inspect the results:
The final verkko-fillet output directory is created next to the input Verkko directory and named {verkko_directory}_verkko_fillet.

```bash
ls verkko_asm_verkko_fillet/figs/
ls verkko_asm_verkko_fillet/stats/
ls verkko_asm_verkko_fillet/chromosome_assignment/
```