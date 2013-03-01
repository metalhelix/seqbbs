context("Basic Validitiy Checks")

test_that("function seqbbs", {
  window <- 20
  threshold <- 0.55
  ratios <- c(0.2, 0.4, 0.6, 0.4)
  test_filename <- system.file("extdata", "test.txt", package="seqbbs")
  ratios <- read.table(test_filename, header = FALSE)
  results <- seqbbs(ratios, window = window, threshold = threshold)
  
})
