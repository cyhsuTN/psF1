psF1: Interval estimation, hypothesis testing, and power and sample size
calculation for F-beta scores and differences in F-beta scores
================
Chih-Yuan Hsu

Sept/10/2025

Chih-Yuan Hsu, Qi Liu, and Yu Shyr (2025). A Unified Framework for Statistical Inference and Power Analysis of Single and Comparative F-beta Scores. Under review.

## Installation

Download psF1_0.4.7.tar.gz and locally install it, or execute the
following code:

``` r
library(devtools)
install_github("cyhsuTN/psF1")
```

## Usage

``` r
library(psF1)
```

## F1 score of a single classifier

Suppose that a test set containing N=58 classification instances and
S=12 positives is used to evaluate the F1 (beta=1) score of a new classifier. The
observed TP=10 (True positive) and FP=1 (False positive).

### Confidence interval

What is the 95% confidence interval for the F1 score of the new
classifier?

``` r
ciFbOne(d=10, b=1, N=58, s=12,
        beta=1,
        alpha=0.05)
```

    ##     Fb.est    95%CI.L    95%CI.U Fb.est.raw 
    ##  0.8400000  0.6486486  1.0000000  0.8695652

### Hypothesis testing

A two-sided test is used to detect whether the F1 score of the new
classifier differs from a target F1 score defined by pp=9/11 (precision, also PPV) and
ps=0.75 (sensitivity, also recall).

``` r
testFbOne(d=10, b=1, N=58, s=12,
          Fb=NULL, pp=9/11, ps=0.75,
          beta=1,
          alternative=c("two.sided", "less", "greater")[1])
```

    ##      fb P-value      Fb      pp      ps 
    ## 0.86957 0.42744 0.78261 0.81818 0.75000

A two-sided test is used to detect whether the F1 score of the new
classifier differs from a target F1 score of 0.783.

``` r
testFbOne(d=10, b=1, N=58, s=12,
          Fb=0.783, pp=NULL, ps=NULL,
          beta=1,
          alternative=c("two.sided", "less", "greater")[1])
```

    ##        fb   P-value        Fb pp.argsup ps.argsup 
    ##   0.86957   0.43398   0.78300   0.83000   0.74104

### Power and sample size calculation

Assume that we plan to re-design such a study.

How large N and S are needed to achieve 80% statistical power to detect
a difference between a new classifier’s F1 score, defined by pp=0.9 and
ps=0.85, and a target F1 score, defined by pp=9/11 and ps=0.75, using a
two-sided test with a 5% significance level?

``` r
powerFbOne(N=round(58*6.5), s=round(12*6.5),
           Fb1=NULL, pp1=0.90, ps1=0.85,
           Fb=NULL, pp=9/11, ps=0.75,
           beta=1,
           alternative=c("two.sided", "less", "greater")[1],
           alpha.level=0.05)
```

    ##   Power      pp      ps      Fb     pp1     ps1     Fb1 
    ## 0.80560 0.81818 0.75000 0.78261 0.90000 0.85000 0.87429

``` r
c(N=round(58*6.5), s=round(12*6.5))
```

    ##   N   s 
    ## 377  78

How large N and S are needed to achieve 80% statistical power to detect
a difference between a new classifier’s F1 score of 0.874 and a target
F1 score defined by pp=9/11 and ps=0.75, using a two-sided test with a
5% significance level?

``` r
powerFbOne(N=round(58*6.8), s=round(12*6.8),
           Fb1=0.874, pp1=NULL, ps1=NULL,
           Fb=NULL, pp=9/11, ps=0.75,
           beta=1,
           alternative=c("two.sided", "less", "greater")[1],
           alpha.level=0.05)
```

    ##   Power      pp      ps      Fb     pp1     ps1     Fb1 
    ## 0.80005 0.81818 0.75000 0.78261 1.00000 0.77620 0.87400

``` r
c(N=round(58*6.8), s=round(12*6.8))
```

    ##   N   s 
    ## 394  82

How large N and S are needed to achieve 80% statistical power to detect
a difference between a new classifier’s F1 score of 0.874 and a target
F1 score of 0.783, using a two-sided test with a 5% significance level?

``` r
powerFbOne(N=round(58*8.2), s=round(12*8.2),
           Fb1=0.874, pp1=NULL, ps1=NULL,
           Fb=0.783, pp=NULL, ps=NULL,
           beta=1,
           alternative=c("two.sided", "less", "greater")[1],
           alpha.level=0.05)
```

    ##   Power      pp      ps      Fb     pp1     ps1     Fb1 
    ## 0.80799 0.95000 0.66594 0.78300 1.00000 0.77620 0.87400

``` r
c(N=round(58*8.2), s=round(12*8.2))
```

    ##   N   s 
    ## 476  98

## Comparison in F1 scores of two classifiers

Suppose that a test set containing N=58 classification instances and
S=12 positives is used to compare the F1 scores of two classifiers. The
observed TP for classifier 1 and 2 are 10 and 9, respectively, and the
observed FP are 1 and 2.

### Confidence interval

What is the 95% confidence interval for the difference between the two
classifiers’ F1 scores?

``` r
ciFbTwo(d1=10, b1=1, 
        d2=9, b2=2,
        N=58, s=12,
        beta=1,
        alpha=0.05)
```

    ##     dFb.est     95%CI.L     95%CI.U dFb.est.raw 
    ##  0.08000000 -0.16879795  0.33333333  0.08695652

### Hypothesis testing

A two-sided test is used to detect whether the difference between the
two classifiers’ F1 scores is significant.

The null distribution is generated using
$F_{1,1}=F_{1,2}=\widehat{F}_{1,0}$

``` r
testFbTwo(d1=10, b1=1,
          d2=9, b2=2,
          Fb="est",
          N=58, s=12,
          beta=1,
          alternative=c("two.sided", "less", "greater")[1])
```

    ##     fb1     fb2     dfb P-value Fb0.est  pp.bar  ps.bar 
    ## 0.86957 0.78261 0.08696 0.47818 0.82609 0.86364 0.79167

The null distribution is generated using $F_{1,1}=F_{1,2}=F_{1,0}$ for
all all possible $F_{1,0}$

``` r
testFbTwo(d1=10, b1=1,
          d2=9, b2=2,
          Fb=NULL,
          N=58, s=12,
          beta=1,
          alternative=c("two.sided", "less", "greater")[1])
```

    ##        fb1        fb2        dfb    P-value Fb0.argsup  pp.argsup  ps.argsup 
    ##    0.86957    0.78261    0.08696    0.67718    0.60000    1.00000    0.42857

### Power and sample size calculation

Assume that we plan to re-design such a comparison.

How large N and S are needed to achieve 80% statistical power to detect
a difference between classifier 1’s F1 score, defined by pp=0.9 and
ps=0.85, and classifier 2’s F1 score, defined by pp=9/11 and ps=0.75,
using a two-sided test with a 5% significance level?

``` r
powerFbTwo(N=round(58*13), s=round(12*13),
           Fb1=NULL, pp1=0.90, ps1=0.85,
           Fb2=NULL, pp2=9/11, ps2=0.75,
           beta=1,
           alternative=c("two.sided", "less", "greater")[1],
           alpha.level=0.05)
```

    ##   Power     pp1     ps1     Fb1     pp2     ps2     Fb2 
    ## 0.80537 0.90000 0.85000 0.87429 0.81818 0.75000 0.78261

``` r
c(N=round(58*13), s=round(12*13))
```

    ##   N   s 
    ## 754 156

