## Running Verkko

Verkko supports multiple assembly configurations, depending on the available sequencing data and phasing resources. For comprehensive documentation, refer to the official Verkko repository: https://github.com/marbl/verkko.

### Basic configuration: HiFi + ONT
Use this configuration when both PacBio HiFi and ONT long reads are available.

```bash
verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz
```

If specific contaminant or auxiliary sequences are known in advance (for example, mitochondrial DNA or rDNA), include `--screen` to produce separate FASTA outputs and representative (exemplar) sequences.

```bash
verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz \
    --screen MT mt.fasta.gz \
    --screen rDNA rDNA.fasta.gz
```

### Trio-informed assembly
If sequencing data from both parents are available (commonly Illumina short reads), parental haplotype-specific k-mers can be used to improve phasing.

```bash
## Build parental hapmer databases (see Verkko documentation for details).

## Run Verkko with trio hapmer support.
verkko -d asm \
  --hifi hifi/*.fastq.gz \
  --nano ont/*.fastq.gz \
  --hap-kmers paternal_compress.k30.hapmer.meryl \
              maternal_compress.k30.hapmer.meryl \
              trio
```

### Hi-C or Pore-C informed assembly
Hi-C or Pore-C data can be incorporated to improve long-range phasing and scaffolding consistency.

```bash
verkko -d asm \
  --hifi hifi/*.fastq.gz \
  --nano ont/*.fastq.gz \
  --hic1 hic/*R1*fastq.gz \
  --hic2 hic/*R2*fastq.gz
```

### Combined Hi-C and trio data
This section will be provided in a future update.

```bash
# To be added
```

### Using ONT-corrected reads
Some workflows (for example, HiFiasm-based preprocessing) can generate ONT-corrected reads. Depending on the experimental design, these reads may be incorporated with HiFi data to increase effective backbone coverage.

```bash
# To be added
```