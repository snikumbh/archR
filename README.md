
Note: _This package is currently under development. So, please bear with me while I put the final blocks together. Thanks for your understanding!_ 

# archR
<!-- badges: start -->
[![Travis build status](https://travis-ci.org/snikumbh/archR.svg?branch=master)](https://travis-ci.org/snikumbh/archR)
<!-- badges: end -->
This is an R package implementing an unsupervised, non-negative matrix factorization-based method (of the same name) for discovery of sequence architectures.


## Installation

### Python scikit-learn dependency
This package requires the Python module scikit-learn. See installation instructions [here](https://scikit-learn.org/stable/install.html).


### To install this package, use 

```r
# install.packages("devtools")
devtools::install_github("snikumbh/archR", build_vignettes = TRUE)
``` 


### Usage
```r
# tssSeqs holds the one-hot encoded input matrix of all sequences
tssSeqs <- prepare_data_from_FASTA(inputFastaFilename)

# Set archR configuration
thisConfig <- archR::archRSetConfig(innerChunkSize = 500,
                                    kMin = 1, kMax = 8, 
				    parallelize = TRUE, nCores = 4,
		                    cvFolds = 3, nIterationsUse = 100)
# Call archR
archRresult <- archR::archR(config = thisConfig, 
                            tss.seqs = tssSeqs)
```


## Troubleshooting Installation

- List some points or link to already known issues reported/resolved
- Possible problems: 
 - installing Python module scikit-learn 
 - getting reticulate to work

# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/archR/issues/new)
