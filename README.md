<!-- badges: start -->
[![Travis build status](https://travis-ci.org/snikumbh/archR.svg?branch=master)](https://travis-ci.org/snikumbh/archR)
<!-- badges: end -->

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/snikumbh/archR/branch/master/graph/badge.svg)](https://codecov.io/gh/snikumbh/archR?branch=master)
<!-- badges: end -->

# archR
R package for discovery of clusters with characteristic sequence architectures using NMF

## Installation

# Python scikit-learn dependency
This package requires the Python module scikit-learn. See installations instructions [here](insert-link).


# To install this package, use 

```
library(devtools)
withr::with_libpaths(<path_here>, 
			install_github("snikumbh/archR", 
				auth_token="29e825452191fa6cdc161719d5a3ddeebe2b017b", 
				build = TRUE)
			)
``` 

If R complains, "Skipping install of 'archR' from a github remote, the SHA1 (7bfe4812) has not changed since last install.
  Use `force = TRUE` to force installation", do as it suggests.

```
withr::with_libpaths(<path_here>, 
			install_github("snikumbh/archR", 
				auth_token="29e825452191fa6cdc161719d5a3ddeebe2b017b", 
				build = TRUE, 
				force = TRUE)
			)
```


## Troubleshooting Installation

- List some points or link to issues if they are already present
- Possible problems: 
 - installing Python module scikit-learn, 
 - getting reticulate to work

# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/archR/issues/new)
