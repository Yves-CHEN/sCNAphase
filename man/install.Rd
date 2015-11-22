\name{How to install}
\alias{sCNAPhase,SCNAPhase}
\docType{package}
\title{
    How to install the tool and address the dependencies
}
\description{
sCNAphase is dependent on the  \href{http://ab-initio.mit.edu/wiki/index.php/NLopt_C-plus-plus_Reference}{NLOPT},
which is a C++ library of local and global optimation algorithms.
The NLOPT library is expected to be installed beforehand.
To locate the header and the library files of NLOPT, the follow two environment variables needs to be set: 

\code{
export PKG_NLOPT_LIBS="/folder/to/NLOPT/lib"    
} \cr
\code{
 export PKG_NLOPT_INCLUDE="/folder/to/NLOPT/include"
}

sCNAphase employs the \href{http://www.boost.org/}{Boost C++ library} for arithmetic calculations.
This dependency can be addressed by installing the BH package from CRAN.
After these two steps, the following command will install the sCNAphase package. \cr
\code{ R CMD INSTALL ./sCNAphase }

For faster speed, sCNAphase enables parallel computing using \href{http://openmp.org/wp/}{OpenMP}, so that it can run on multiple CPUs.
This installation procedure has been tested on Linux with GCC-4.4.3.

}

\note{This compiling procedure was tested on centos, ubuntu with GCC-4.4.3, R-3.0.1.}
