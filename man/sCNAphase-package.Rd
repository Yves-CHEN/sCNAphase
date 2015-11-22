\name{Package overview}
\alias{sCNAPhase,SCNAPhase, sCNAphase, sCNAphase-package}
\docType{package}
\title{
    sCNAphase: somatic copy number profiling based on allelic phases. 
}
\description{
sCNAphase is an R package, designed for estimating haploid somatic copy number alternations (sCNAs), tumor cellularity and tumor ploidy, 
based on whole genome sequencing (WGS) or whole exome sequencing (WES) data. To detect the somatic alterations and avoid the GC bias,
a patient-matched normal samples is required and expected to be sequenced based on the same protocal and platform as the tumor samples. 

To generate the input for sCNAphase, a pre-processing step is expected to produce two sets of vcf files and a set of haps files.
 \enumerate{
    \item Set N:  a vcf file with germ-line SNPs called from a normal sample for each chromosome.  
    \item Set T:  a vcf file called from a tumor sample at the germ-line SNPs for each chromosome. 
    \item Set H:  a hap file with the phase information for each germ-line SNP for each chromosome.   
  }
Set N and Set T is generated from samtools mpileup.   \cr
Set C is calculated based on Set A using \href{https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html}{SHAPEIT}.

Given these files, sCNAphase generates the sCNAs profile at haploid level, the tumor purity and ploidy estimations through the inferCNA function. The sCNAs profile can be formated into segmentation files, vcf files or visualized in d.SKY plots.

}
\author{
Wenhan CHEN
}
\references{
Wenhan CHEN, sCNAphase package, \url{https://github.com/Yves-CHEN/sCNAphase}
}
\keyword{ copy number }
\keyword{ somatic }
\keyword{ tumor purity }
\keyword{ tumor cellularity }
\keyword{ tumor ploidy }
\keyword{ haploid copy number }
\keyword{ allelic phases }
\seealso{
\code{\link{inferCNA}}
\code{\link{genSegFile}}
\code{\link{produceDSKY}}
}

\note{
Before running the following code in command line, specify the the number of CPUs by \cr
\code{ export OMP_NUM_THREADS=2   }\cr
so that this would run on 2 CPUs.
}

\examples{
# -----------------------------------------------------------------------------
# This is a demo written in R, with minimal number of paramenters.
# -----------------------------------------------------------------------------
anaName = "inferCN"  # This specifies the name of the analysis.
    # It can be any name, but better to be meaningful and unique,
    # as the result files are made up of anaName.

chroms  = c(1:22)    # This specify chromosomes of genome to consider.  
    #1:22 means the 22 autosomes.

nPrefix = "/baseDir1/filename" # inferCNA function will try to locate
    # the a vcf for a normal genome, 
    # named as /baseDir1/filename.chr{1...22}.vcf,   (Set N)
    # a hap file /baseDir1/filename.chr{1...22}.haps (Set H)

tPrefix = "/baseDir2/filename"  # inferCNA will try to locate the a vcf
    #for a tumor genome, named as /baseDir2/filename.chr{1...22}.vcf (Set T)

inferCNA (anaName, nPrefix, tPrefix, chroms)   # inferCNA will profile 
    # sCNAs and generate a R dat file called res.{anaName}.phased.chr.W.dat
    # in the current  directory, based on which sCNAphase can then format
    # the estimation to the segmentation file, the d.SKY plot and a vcf file.
                                             
genSegFile(anaList= anaName, outdir = "test")   # This generates a *.csv file.
    # Each row corresponds to a particular sCNAs with chrID, start, end, 
    # copy number, allelic copy number. 

                                                
produceDSKY(anaList= anaName, outDir = "test")   
# This will generate the d.SKY plot into a pdf file.                                               

}


