context("Plotting features matrix as ggheatmap")
library(ggplot2)
library(reshape2)
library(dplyr)


test_that("Matrix has 4 rows", {

  testPwmMat <- matrix(rnorm(100), nrow = 2)
  testPositionLabels <- seq(25)
  expect_error(plot_ggheatmap(testPwmMat, position_labels = testPositionLabels)
               , "Expecting a matrix with 4 rows corresponding to DNA alphabet")

})

test_that("Given object is matrix", {

  testPwmMat <- rnorm(200) # err
  testPositionLabels <- seq(25)
  expect_error(plot_ggheatmap(testPwmMat, position_labels = testPositionLabels)
               , "Expecting a matrix with 4 rows")

})

test_that("Position labels inadequate", {

  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(20)
  expect_error(plot_ggheatmap(testPwmMat, position_labels = testPositionLabels)
               , "Inadequate")

})

test_that("Position labels over-abundant", {

  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(50)
  expect_error(plot_ggheatmap(testPwmMat, position_labels = testPositionLabels)
               , "Overabundant")

})

# test_that("ggheatmap plotting works", {
#
#   # test
#   testPositionLabels <- seq(5)
#   testPwmMat <- matrix(rnorm(20), nrow=4)
#   rownames(testPwmMat) <- c('a', 'c', 'g', 't')
#   pwmMat_df <- as.data.frame(testPwmMat)
#   #
#   pwmMat_df <- mutate(pwmMat_df, Nucleotides = rownames(pwmMat_df))
#   colnames(pwmMat_df) <- c(testPositionLabels, "Nucleotides")
#   #
#   pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars=c("Nucleotides"), variable.name = "positions")
#
#   p1 <- ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(x = positions,
#                                                              y = Nucleotides,
#                                                              fill = value)) +
#                 geom_tile() +
#                 theme_bw() +
#                 xlab(label = "Positions") +
#                 scale_fill_gradient2(name = "", #"Loading",
#                                      low = "white", #"#FFFFFF",
#                                      mid = "white", #FFFFFF",
#                                      high = "#012345"
#                 ) +
#                 theme(legend.position = "top",
#                       legend.justification = "center"
#                 )
#   #
#   expect_equal(plot_ggheatmap(testPwmMat, position_labels = testPositionLabels), p1)
#
# })
