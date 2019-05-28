context("Plotting seqlogo w/ ggseqlogo")
library(ggseqlogo)


test_that("Matrix has 4 rows", {

    testPwmMat <- matrix(rnorm(100), nrow = 2)
    testPositionLabels <- seq(25)
    expect_error(plot_ggseqlogo(testPwmMat, position_labels = testPositionLabels)
                 , "Expecting a matrix with 4 rows corresponding to DNA alphabet")

})

test_that("Given object is matrix", {

  testPwmMat <- rnorm(200) # err
  testPositionLabels <- seq(25)
  expect_error(plot_ggseqlogo(testPwmMat, position_labels = testPositionLabels)
               , "Expecting a matrix with 4 rows")

})

test_that("Position labels inadequate", {

  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(20)
  expect_error(plot_ggseqlogo(testPwmMat, position_labels = testPositionLabels)
               , "Inadequate")

})

test_that("Position labels over-abundant", {

  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(50)
  expect_error(plot_ggseqlogo(testPwmMat, position_labels = testPositionLabels)
               , "Overabundant")

})


