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
heterozygous deletions. An example of this can be found at (url N.A.). This allows visual inspections
of chromsomal, focal copy number deletions, gains or amplification. Generation of the cytobands requires 
the quantsmooth R packages.
}

\usage{ genSegFile(anaList= anaName, outdir = "test") 
produceDSKY(anaList= anaName, outdir = "test") }
\arguments{
\item{anaName}{Any label for the analysis. This should be the same as the anaName specified for function inferCNA}
\item{outdir}{The output directory.}

}
\examples{
genSegFile(anaList= anaName, outdir = "test")   # This generates a *.csv file.
    # Each row corresponds to a particular sCNAs with chrID, start, end, 
    # copy number, allelic copy number. 
                                                
produceDSKY(anaList= anaName, outDir = "test")   
# This will generate the d.SKY plot into a pdf file.                                               
}
