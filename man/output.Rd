\name{Output Formating}
\alias{genSegFile, produceDSKY}
\title{Generate CN segmentation, vcf files and d.SKY plot
}
\alias{inferCNA}

\description{
By default, sCNAphase generates a R data file *.dat which contains all the information.
This raw result can be formatted to more readable outputs using genSegFile, produceDSKY.

The genSegFile produces a *.csv file, similar to *.bed file format, which includes 5 columns:
    

  \tabular{rlllll}{
    \tab chr \tab    start  \tab    end \tab   CN \tab    mCN
  }

which stands for chromsome identifier, start position, end position, the overall copy number, the copy number for the less amplified allele respectively. If cases of LOHs and heterozygous deletions, mCN is 0.

The digital SKY plot looks like spectural karyotypiing (SKY) plots, which shows the chromosomes 
in cytobands and marks the copy number for each regions, the regions with LOHs, the regions with
heterozygous deletions. An example of this can be found at (https://figshare.com/authors/_/1365237). This allows visual inspections
of chromsomal, focal copy number deletions, gains or amplification. Generation of the cytobands requires 
the quantsmooth R packages.
}

\usage{ genSegFile(anaList= anaName, outdir = "test",
                     mlRemoveLevel = 0.04, ifload = T) 
produceDSKY(anaList= anaName, outdir = "test",
            mlRemoveLevel = 0.04, ifload = T) }
\arguments{
\item{anaName}{Any label for the analysis. This should be the same as the anaName specified for function inferCNA}
\item{outdir}{The output directory.}
\item{mlRemoveLevel}{When merging phased allelic depths into PHF, the merging error can rise at the boundary between two neighbor sCNAs. mlRemoveLevel determines to which degree, the PHFs needs to be filtered due to merging error. By default, the value is 0.04, so that 4\% of the PHFs will be removed. }
\item{ifload}{By default, the anaName allows the function to identify the location of a *.dat. Then within the functions, the *.dat file is loaded. However, if the ifload is set to F, the *.dat file can be specified and loaded from outside this function. }

}
\examples{
genSegFile(anaList= anaName, outdir = "test")   # This generates a *.csv file.
    # Each row corresponds to a particular sCNAs with chrID, start, end, 
    # copy number, allelic copy number. 
                                                
produceDSKY(anaList= anaName, outDir = "test")   
# This will generate the d.SKY plot into a pdf file.                                               
}
