

# archR
<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![DOI](https://zenodo.org/badge/188449833.svg)](https://zenodo.org/badge/latestdoi/188449833)
[![Build status](https://travis-ci.org/snikumbh/archR.svg?branch=master)](https://travis-ci.org/snikumbh/archR)
[![Codecov test coverage](https://codecov.io/gh/snikumbh/archR/branch/master/graph/badge.svg)](https://codecov.io/gh/snikumbh/archR?branch=master)
[![R build status](https://github.com/snikumbh/archR/workflows/R-CMD-check/badge.svg)](https://github.com/snikumbh/archR/actions)
<!-- badges: end -->

Note: _This package is currently under development. So, please bear with me while I put the final blocks together. Thanks for your understanding!_ 


archR is an unsupervised, non-negative matrix factorization (NMF)-based algorithm for discovery of sequence architectures de novo.
Below is a schematic of archR's algorithm.

<img src="https://github.com/snikumbh/archR/blob/master/vignettes/archR_algorithm_1080p_cropped.gif" width="550" align="center">


## Installation

### Python scikit-learn dependency
This package requires the Python module scikit-learn. Please see installation instructions [here](https://scikit-learn.org/stable/install.html).


### To install this package, use 

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")   
}

remotes::install_github("snikumbh/archR", build_vignettes = FALSE)
``` 



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
                                               raw_seq = TRUE)

nSeqs <- length(inputSeqsRaw)
positions <- seq(1, Biostrings::width(inputSeqsRaw[1]))

# Set archR configuration
# Most arguments have default values
archRconfig <- archR::archR_set_config(
        parallelize = TRUE,
        n_cores = 2,
        n_runs = 100,
        k_min = 1,
        k_max = 20,
        mod_sel_type = "stability",
        bound = 10^-6,
        chunk_size = 100,
	result_aggl = "ward.D",
	result_dist = "euclid",
        flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
                     plot = FALSE)
        )

#
### Call/Run archR
archRresult <- archR::archR(config = archRconfig,
                               seqs_ohe_mat = inputSeqsMat,
                               seqs_raw = inputSeqsRaw,
                               seqs_pos = positions,
                               total_itr = 2,
			       set_ocollation = c(TRUE, FALSE))

```


# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/archR/issues/new)
