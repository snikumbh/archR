

# archR
<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Build status](https://travis-ci.org/snikumbh/archR.svg?branch=master)](https://travis-ci.org/snikumbh/archR)
[![Codecov test coverage](https://codecov.io/gh/snikumbh/archR/branch/master/graph/badge.svg)](https://codecov.io/gh/snikumbh/archR?branch=master)
<!-- badges: end -->

Note: _This package is currently under development. So, please bear with me while I put the final blocks together. Thanks for your understanding!_ 


archR is an unsupervised, non-negative matrix factorization (NMF)-based algorithm for discovery of sequence architectures de novo.
Below is a schematic of archR's algorithm.


<img src="reference/figures/archR_algorithm_1080p_cropped.gif" width="550" align="center">


## Installation

### Python scikit-learn dependency
This package requires the Python module scikit-learn. See installation instructions [here](https://scikit-learn.org/stable/install.html).


### To install this package, use 

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")   
}

remotes::install_github("snikumbh/archR", build_vignettes = FALSE)
``` 

If the above command produces errors, see this troubleshooting [note](https://github.com/snikumbh/archR#troubleshooting-installation) below.


### Usage
```r
# load package
library(archR)
library(Biostrings)


# Creation of one-hot encoded data matrix from FASTA file
# You can use your own FASTA file instead
inputFastaFilename <- system.file("extdata", "example_data.fa", 
                                  package = "archR", 
                                  mustWork = TRUE)

# Specifying dinuc generates dinucleotide features
inputSeqsMat <- archR::prepare_data_from_FASTA(inputFastaFilename,
                                                  sinuc_or_dinuc = "dinuc")

inputSeqsRaw <- archR::prepare_data_from_FASTA(inputFastaFilename, 
                                               rawSeq = TRUE)

nSeqs <- length(inputSeqsRaw)
positions <- seq(1, Biostrings::width(inputSeqsRaw[1]))
sinuc <- Biostrings::DNA_BASES

# Set archR configuration
# Most arguments have default values
archRconfig <- archR::archRSetConfig(
        parallelize = TRUE,
        nCoresUse = 4,
        nIterationsUse = 100,
        kMin = 1,
        kMax = 20,
        modSelType = "stability",
        tol = 10^-4,
        bound = 10^-8,
        innerChunkSize = 500,
        flags = list(debugFlag = FALSE, timeFlag = TRUE, verboseFlag = TRUE,
                     plotVerboseFlag = FALSE)

#
### Call/Run archR
perform_iters <- 2
archRresult <- archR::archR(config = archRconfig,
                               seqsMat = inputSeqsMat,
                               seqsRaw = inputSeqsRaw,
                               seqsPositions = positions,
                               thresholdItr = perform_iters)

```

## Troubleshooting Installation

1. Error during installation (with `remotes::install_github` or `devtools::install_github`. If the error occurs during `** testing if installed package can be loaded`, try adding the `--no-test-load` option as `remotes::install_github("snikumbh/archR", INSTALL_opts = c("--no-test-load"))`


# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/archR/issues/new)
