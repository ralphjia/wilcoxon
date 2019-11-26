# wilcoxon

 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/ralphjia/wilcoxon.svg?branch=master)](https://travis-ci.org/ralphjia/wilcoxon)
  <!-- badges: end -->

<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/ralphjia/wilcoxon/branch/master/graph/badge.svg)](https://codecov.io/gh/ralphjia/wilcoxon?branch=master)
  <!-- badges: end -->

This package contains a function wilcoxon_test(), which can perform a Wilcoxon signed-rank test or a Wilcoxon rank sum test, also known as a Mann-Whitney test. 

The Wilcoxon signed rank test is a non-parametric hypothesis test that can be applied to either a single sample or a pair of samples of equal size. The Wilcoxon signed rank test can be used to test the null hypothesis that the median of a distribution is equal to some value mu, or that the median of the difference of two distributions is equal to mu. This test can be used in place of a one-sample t-test or a paired t-test when the data cannot be assumed to be normally distributed.

The Wilcoxon rank sum test is another non-parametric hypothesis test which can be applied to two samples, of any size. The Wilcoxon rank sum test can be used to test the null hypothesis that the difference of medians is equal to some location shift, mu. For this test, the two samples must be independent of each other, but it also does not require that the data be normally distributed. 

Both of these tests are built into R through the function wilcox.test(). A comparison of wilcoxon_test() vs. wilcox.test() with regards to accuracy and efficiency can be found in the vignetttes. wilcox.test() is a more sophisticated function, in that it allows for computation of confidence intervals, whereas wilcoxon_test() will only calculate a test statistic and p-value.

# Installation
```{r}
devtools::install_github("ralphjia/wilcoxon")
```

