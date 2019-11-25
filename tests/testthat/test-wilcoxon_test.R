context("test-wilcoxon_test")

test_that("One sample Signed Rank test works", {
  expect_equal(wilcoxon_test(1:5), wilcox.test(1:5))
  expect_equal(wilcoxon_test(-5:-1), wilcox.test(-5:-1))
  expect_equal(wilcoxon_test(3:7, mu = 4), wilcox.test(3:7, mu = 4, exact = FALSE))
  expect_equal(wilcoxon_test(1:5, alternative = "l"), wilcox.test(1:5, alternative = "l"))
  expect_equal(wilcoxon_test(1:5, alternative = "g"), wilcox.test(1:5, alternative = "g"))
  expect_equal(wilcoxon_test(-45:55), wilcox.test(-45:55))
  expect_equal(wilcoxon_test(-45:55, alternative = "l"), wilcox.test(-45:55, alternative = "l"))
  expect_equal(wilcoxon_test(-45:55, alternative = "g"), wilcox.test(-45:55, alternative = "g"))
})

test_that("Paired Signed Rank test works", {
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T), wilcox.test(1:5, 3.5:-0.5, paired = T))
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T, alternative = "l"), wilcox.test(1:5, 3.5:-0.5, paired = T, alternative = "l"))
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T, alternative = "g"), wilcox.test(1:5, 3.5:-0.5, paired = T, alternative = "g"))
})

test_that("Paired test error if y is missing", {
  expect_error(wilcoxon_test(1:5, paired = T), "'y' is missing for paired test")
})

test_that("Paired test error if x and y have different lengths", {
  expect_error(wilcoxon_test(1:5, 1:4, paired = T), "'x' and 'y' must have the same length")
})

test_that("Warning if exact is true but there are ties", {
  expect_warning(wilcoxon_test(-2.5:2.5, exact = T), "cannot compute exact p-value with ties")
})

test_that("Warning if exact is true but there are zeroes", {
  expect_warning(wilcoxon_test(0:5, exact = T), "cannot compute exact p-value with zeroes")
})

test_that("Wilcoxon Rank-Sum test works", {
  expect_equal(wilcoxon_test(1:5, 0.5:8.5), wilcox.test(1:5, 0.5:8.5))
  expect_equal(wilcoxon_test(1:5, 0.5:8.5, alternative = "l"), wilcox.test(1:5, 0.5:8.5, alternative = "l"))
  expect_equal(wilcoxon_test(1:5, 0.5:8.5, alternative = "g"), wilcox.test(1:5, 0.5:8.5, alternative = "g"))
  expect_equal(wilcoxon_test(0.5:8.5, 1:5), wilcox.test(0.5:8.5, 1:5))
  expect_equal(wilcoxon_test(-20:30, -25.5:20.5), wilcox.test(-20:30, -25.5:20.5))
  expect_equal(wilcoxon_test(-20:30, -25.5:20.5, alternative = "l"), wilcox.test(-20:30, -25.5:20.5, alternative = "l"))
  expect_equal(wilcoxon_test(-20:30, -25.5:20.5, alternative = "g"), wilcox.test(-20:30, -25.5:20.5, alternative = "g"))
})

test_that("Warning if exact is true but there are ties (for rank-sum test)", {
  expect_warning(wilcoxon_test(1:5, 4:7, exact = T), "cannot compute exact p-value with ties")
})


