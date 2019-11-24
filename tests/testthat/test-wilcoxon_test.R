context("test-wilcoxon_test")

test_that("multiplication works", {
  expect_equal(wilcoxon_test(1:5)$statistic, wilcox.test(1:5)$statistic)
})
