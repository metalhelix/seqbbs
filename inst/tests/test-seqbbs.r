context("Basic Validitiy Checks")

test_that("function seqbbs", {
  window <- 20
  threshold <- 0.70
  #ratios <- c(0.2, 0.4, 0.6, 0.4)
  test_filename <- system.file("extdata", "paper.txt", package="seqbbs")
  ratios <- read.table(test_filename, header = FALSE)
  results <- seqbbs(ratios, window = window, threshold = threshold)
  thresholded_changepoints <- changepoints(results, threshold = threshold)
  
  expected_results_filename <- system.file("extdata", "expected_results.csv", package="seqbbs")
  expected <- read.table(expected_results_filename, header = TRUE, sep = ",")
  
  expect_that(thresholded_changepoints$loci, equals(expected$loci))
  
  ratios <- round(thresholded_changepoints$mean_ratio, digits = 4)
  expect_that(ratios, is_equivalent_to(expected$mean_ratio))
  
  low <- round(thresholded_changepoints$confidence_lower_bound, digits = 4)
  expect_that(low, is_equivalent_to(expected$confidence_lower_bound))
  
  hi <- round(thresholded_changepoints$confidence_upper_bound, digits = 4)
  expect_that(hi, is_equivalent_to(expected$confidence_upper_bound))
  
})
