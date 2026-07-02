## Check inputs for Verkko


### Coverage check
Verkko is designed to use both PacBio HiFi and ONT long-read data. Before assembly, verify that each dataset provides sufficient depth. As a practical baseline, we recommend at least **40x** coverage for both HiFi and ONT, with **60x or higher** preferred for improved assembly continuity and phasing performance.

Coverage can be estimated as:

$$
	ext{coverage} = \frac{\text{total sequenced bases}}{\text{haploid genome size (bp)}}
$$

where total sequenced bases is the sum of read lengths across all reads in the dataset.

You can use script below: 
```bash
# to be added
```