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

test_that("Handling empty matrix", {

  testPwmMat <- matrix()
  testPositionLabels <- seq(25)
  expect_error(plot_ggseqlogo(testPwmMat, position_labels = testPositionLabels)
               , "Empty")

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

test_that("ggseqlogo plotting works", {
  # setting seed enables proper comparison between ggplot objects since we use
  # rnorm
  set.seed(11223344)
  # test variables
  testPositionLabels <- seq(25)
  testPwmMat <- matrix(rnorm(100), nrow=4)
  testPwmMat <- make_sinuc_PWMs(testPwmMat)
  p1 <- plot_ggseqlogo(pwmMat = testPwmMat, position_labels =  testPositionLabels)
  # test plot
  vdiffr::expect_doppelganger("ggseqlogo plot example", p1)

})
