library(testthat)
library(ggplot2)

# Assemble
mix <- c(0.1, 0.2, 0.3, 0.4)
depths <- c(10, 20, 30, 40)
df_barcodes <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1), nrow = 4, byrow = TRUE)
muts <- c("mut1", "mut2", "mut3")
Y <- c(0.1, 0.2, 0.3, 0.4)
lmps <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1), nrow = 4, byrow = TRUE)

# Tests for Freyja model
test_that("freyja model produces non-zero coefficients", {
  # Act
  freyja_coeffs <- provoc::freyja(mix, depths, df_barcodes, muts)

  # Assert
  expect_true(any(freyja_coeffs > 0), info = "Freyja model should produce some non-zero coefficients")
})

# Tests for Alcov model
test_that("alcov model produces expected coefficients", {
  # Act
  alcov_coeffs <- provoc::alcov(Y, lmps, muts)

  # Assert
  expected_alcov_coeffs <- c(0.0, 0.3333333, 0.6666667)
  expect_equal(alcov_coeffs, expected_alcov_coeffs, tolerance = 1e-5, 
               info = "Alcov model coefficients should match expected values")
})

# Plotting coefficients for visual comparison
test_that("plot coefficients for visual comparison", {

  freyja_coeffs <- provoc::freyja(mix, depths, df_barcodes, muts)
  alcov_coeffs <- provoc::alcov(Y, lmps, muts)

  coefficients_df <- data.frame(
    model = rep(c("Freyja", "Alcov"), each = 3),
    mutation = rep(muts, 2),
    coefficient = c(freyja_coeffs, alcov_coeffs)
  )
  
  p <- ggplot(coefficients_df, aes(x = mutation, y = coefficient, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Coefficient Comparison between Freyja and Alcov Models",
         y = "Coefficient",
         x = "Mutation") +
    scale_fill_manual(values = c("Freyja" = "blue", "Alcov" = "red")) +
    theme_minimal()
  
  print(p)
})

