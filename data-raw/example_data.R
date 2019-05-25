## code to prepare `example_data` dataset goes here

require(seqinr)
inputFastaFilename <- file.path('data-raw/samarth-dataset-2-4clusters-ZeroMutation.fa')
example_data <- archeR::prepare_data_from_FASTA(inputFastaFilename)
usethis::use_data(example_data)
