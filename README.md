# GenePy_LRSEQ
A gene-based pathogenicity score based on phased variants from long read sequencing

Prerequisites:

  SLURM: Ensure you have SLURM installed and configured on your system.
  
  bcftools: A set of utilities for variant calling and manipulating VCF and BCF files.
  
  conda: A package and environment management system.
  
  CADD: Combined Annotation Dependent Depletion, a tool for scoring the deleteriousness of single nucleotide variants.
  
  VEP: Variant Effect Predictor, a tool for annotating and predicting the effects of variants on genes.
  
  Python: Ensure Python is installed.
  
  NumPy: A library for numerical computations.
  
  Pandas: A library for data manipulation and analysis.
  
  Numba: A library for JIT compiling to optimize numerical functions.
  
  PyArrow: A library for reading and writing data in columnar format.
  
  CUDA: Required if using GPU for computation.

Main scripts:

  pre_1.sh, is designed to process the vcf file from LRSeq data as phased variants are represented with a phase set information
  
  vcf2meta.sh, is designed to process VCF files and generate metadata files with annoation of functional features as input of the genepy algorithm
  
  score.sh, is to process the meta information per gene and generate the pathogenecity score for each haplotype of the gene
  

  Usage:
  
  Navigate to the working directory, sbatch pre_1.sh $input.vcf
  
  sbatch vcf2meta.sh
  
  split gene.list -d -l 800
  
  ls -1 x* | while read i; do sbatch score.sh $i $CADD_CUTOFF;done



  


