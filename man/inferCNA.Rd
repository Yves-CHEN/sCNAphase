\name{Infer copy number}
\alias{inferCNA}
\title{
    A function for haplotype-based allelic copy number alteration inference.
}
\alias{inferCNA}

\description{This function performs haplotype-based allelic copy number alteration inference on the paired normal-tumor sequencing data.
The preprocessing is needed to call the variances and produce the haplotypes.
This function has no return values, but instead produces a *.dat file containing all the information
about the estimation. This *.dat file can be further processed into a vcf files,
segmentation file or a d.SKY plot.


The vcf files and the file with the phase information has to be named in the format of $tPrefix.chr*.vcf and $Prefix.chr*.haps and placed in the same folder.
}



\usage{
inferCNA  <- function(anaName, nPrefix, tPrefix, chroms, doPhase = T, forceRead=F, maxCopyNum=12, mlen = 30, maxiter = 1, ploidy = seq(1, 2.5, 0.1), allelicMapability = F, generateLog = T)
}
\arguments{

   \item{anaName}{Any label for the analysis. The anaName appears in the name of output files.}
   \item{nPrefix}{The prefix of the normal sample files in vcf format.
                     The prefix can include the path specification to the vcf files.
                     For example, nPrefix = "../Data/HCC1143.normal"}
   \item{tPrefix}{The prefix of the tumor sample files in vcf format. Same as above.}
   \item{chroms}{The chromosomes to be included in analaysi.
                     By defaut, chroms=1:22 which represents all the autosomes.}
   \item{doPhase}{If this value is TRUE, the phase information will be included.  Default value is T.
                       The analysis will be based on haplotypes.}
   \item{forceRead}{If this value is F,
                       the function will load the depths information from a temporary *.dat file.
                       If the *.dat doesn't exist, the function will read the depth information
                       from the vcf files specified by nPrefix and tPrefix and generate a temporary *.dat file.
                       Default value is F, since parsing vcf file can cost time. 
                       If it is T, the function will disregard these temporary *.dat file}
   \item{maxCopyNum}{The upper limit of the somatic copy number alteration from normal.
                       Default value is 12. The lower limit is 0 for genotype of homozygous deletions.}
   \item{mlen}{By default, the sCNAphase performs a merge on the SNPs.
                       This specifies how many allelic depths will be merged.
                       Default value is 30 for every 30 SNPs.
                       The greater improves the power of the model for estimation, but reduces the resolution.
                    }
   \item{maxiter}{The function performs an EM estimation.
                       This value specifies the number of iterations. Default value is 1.}
   \item{ploidy}{This specifies a list of possible ploidy values to search,
                       since the average ploidy for cancer is often unknown. sCNAphase chooses the value that
                       maximize the likelihoood function. The value of 1 corresponds to ploidy index of 1, which is
                       ploidy of 2. The default values is 1 to 2.5 incremented by 0.1.  
                       However when a rough value of the average is available, this prior knowledge can be fed to 
                       sCNAphase through this parameter.
                    }
   \item{allelicMapability}{ Allelic bias can cause the imbalance of the allelic depth, even when the 
                                  copy number of reference allele equals to the alternative allele. 
                                  When the allelicMapability is T, sCNAphase will adjust this bias
                                  by correcting the number of reads mapped to each allele.
                      }

   \item{generateLog}{ If this value is TRUE, sCNAphase will generate a pdf report about the analysis. 
                        The default value is T.
                      }


}

\note{Many of the parameters provided can just keep the default value. The anaName, nPrefix, tPrefix are
neccessary, and required to be specified each time. The ploidy parameter is useful, when a rough value of the average ploidy is known.}
