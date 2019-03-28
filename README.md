# Identify_lowVAF_bulk_muts

Identify mutations that are clonally present in a clonally outgrown cell, but subclonally present in a bulk sample. This is done by comparing clones with each other. If a variant is present in a subset of clones, then it must have occurred as a mutation in an ancestor of the clones. Mutations identified in this way can be used to generate phylogenetic trees.

## USAGE:
This script takes a multi-sample vcf as its input. This vcf should contain a bulk sample and at least two clonal samples. Before using this script the vcf should first be filtered with our somaticFilter.py script, or an equivalent filtering method. This script assumes normal karyotypes are present.

```bash
Rscript Determine_lowVAF_bulk.R --vcf Path/To/Filter_somatic_output.vcf --bulk bulk_name --sample_name sample-name --gender [M|F] --out_dir out_dir

```
## Dependencies
- R >= 3.5.1

#### R packages
- tidyverse >= 1.2.1
- VariantAnnotation >= 1.26.1
- GenomicRanges >= 1.32.7
- optparse >= 1.6.1
