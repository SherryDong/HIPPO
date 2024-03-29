# HIPPO V1.0

HIPPO (**H**aplotype **I**m**P**utation for **P**seudogene-mediated c**O**nversion events based on short-read next-generation sequencing data) is a solution started from the general short-read NGS BAM files and dedicated to imputing the genomic haplotypes of functional/pseudogene pairs that have highly-homologous sequences and gene conversion events frequently occurred, along with result visualization including imputed haplotype, the proportion of gene recombination events between two neighboring informative mutations, and specific reads information. 


# Require

```seqtk >= 1.3-r106```

```bedtools >= v2.27.1```

```muscle >= v3.8.31```

```samtools >= v1.8```

```blat >= v.36```

```R >= 3.6.0```

# Demo

HIPPO performs haplotype imputation followed by three main steps: (1) prepare gene-specific reference, (2) generate reads-region mapping content matrix, and (3) functional/pseudogene haplotype imputation and result visualization 

### Demos script for all of the process.
demo.sh

### Demos script for each of the three steps (each script could run separately). 

+ Part I: prepare gene-specific reference

  Step1_demo.sh

+ Part II: generate reads-region mapping content matrix

  Step2_demo.sh
  
  # 

+ Part III: functional/pseudogene haplotype imputation and result visualization 

  Step3_demo.R
  
  Rscript src/pipeline_draw.R
  
  usage: Rscript src/pipeline_draw.R <read-region matrix> <combine_mod_region.bed> <combine_mod.info> <interested position> <output pdf file name> <output txt file name>

# Reference

