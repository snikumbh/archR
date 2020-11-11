# Changes in version 0.1.2 (2020-11-11)

## Enhancements
  - Fix staged installation error. This enables seamless installation on macOS
  - CI with GHA
  - archR result list has new elements:
    - `timeInfo` if timeFlag is TRUE;
    - `clustSol` storing combined clusters from the last iteration of archR.
  - Architectures/clusters sequence logos can be plotted with `auto` y-axis limits
  for information content (0-max instead of 0-2).
  - TFBSTools moved from dependency to suggests. This eases installation of 
  archR by reducing possibility of failed installation.


# Changes in version 0.99.3 (2019-10-21)

## Enhancements
  - Accepts FASTA sequences as DNAStringSet object(s)

