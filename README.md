#  sCNAphase
sCNAphase is an R package, designed for esitmating tumor copy number profile, tumor cellularity, tumor ploidy, based on whole genome sequencing (WGS) or whole exome sequencing (WES) data.
# Requirements
  * R >= (3.1.2)
  * Rcpp
  * BH
  * NLOPT
  


# Install
R CMD INSTALL ./sCNAphase

This only works when the Rcpp, BH, NLOPT are installed. Find the user manual for more information.

# Prepare data
To run the following R code, sCNAphase expects two sets of vcf files and a set of haps files.
   * Set A: a vcf file with germ-line SNPs called from a normal sample for each chromosome.
   * Set B: a vcf file called from a tumor sample at the germ-line SNPs for each chromosome.
   * Set C: a hap file with the phase information for each germ-line SNP for each chromosome.   


> Set A and Set B is generated from samtools mpileup.  
> Set C is calculated based on Set A using [SHAPEIT] (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html).

# Usage
```R
# This is a demo written in R, with minimal number of paramenters.
# -----------------------------------------
anaName = "inferCN"            # Any name. This will make up a part of the result file names.
chroms  = c(1:22)              # The chromosomes that makes up the genome. 1:22 means the 22 autosomes.
nPrefix = "/baseDir/filename"  # inferCNA will try to locate the a vcf from a normal genome,
                               #    named as /baseDir/filename.chr{1...22}.vcf, 
                               #    a hap file /baseDir/filename.chr{1...22}.haps
tPrefix = "/baseDir/filename"  # inferCNA will try to locate the a vcf from a normal genome,
                               #    named as /baseDir/filename.chr{1...22}.vcf
inferCNA (anaName, nPrefix, tPrefix, chroms)  # This will generate a R dat file in the current 
                                              #   directory  called res.{anaName}.phased.chr.W.dat,
                                              #   based on which sCNAphase can then format the
                                              #   estimation to the segmentation file, the d.SKY
                                              #   plot and a vcf file.

genSegFile(anaList= anaName, outdir = "test")   # This generates a *.csv file. Each row corresponds
                                                #   to a particular sCNAs with chrID, start, end, 
                                                #   copy number, allelic copy number. 
produceDSKY(anaList= anaName, outDir = "test")  # The will generate the d.SKY plot into a pdf file.


```
# License
sCNAphase is licensed under L-GPL.
