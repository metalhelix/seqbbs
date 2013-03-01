context("Basic Validitiy Checks")

test_that("function seqbbs", {
  window <- 2
  threshold <- 0.55
  ratios <- c(0.2, 0.4, 0.6, 0.4)
  results <- seqbbs(ratios, window = window, threshold = threshold)
  
})
