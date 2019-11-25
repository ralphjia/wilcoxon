context("test-wilcoxon_test")

test_that("One sample Signed Rank test works", {
  expect_equal(wilcoxon_test(1:5), wilcox.test(1:5))
  expect_equal(wilcoxon_test(1:5, alternative = "l"), wilcox.test(1:5, alternative = "l"))
  expect_equal(wilcoxon_test(1:5, alternative = "g"), wilcox.test(1:5, alternative = "g"))
})

test_that("Paired Signed Rank test works", {
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T), wilcox.test(1:5, 3.5:-0.5, paired = T))
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T, alternative = "l"), wilcox.test(1:5, 3.5:-0.5, paired = T, alternative = "l"))
  expect_equal(wilcoxon_test(1:5, 3.5:-0.5, paired = T, alternative = "g"), wilcox.test(1:5, 3.5:-0.5, paired = T, alternative = "g"))
})
