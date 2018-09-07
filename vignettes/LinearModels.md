---
title: "Crash course in linear modelling for single and multiple genes"
author: "Malte Thodberg"
date: "2018-09-07"
output:
  ioslides_presentation:
    smaller: true
    highlight: tango
    transition: faster
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{EDA}
  %\usepackage[UTF-8]{inputenc}
---

## Introduction

This presentations cover some basic intutions about linear models of gene expression.

We will go through how `lm()` works, focusing on the specification of design or model matrices for some common designs for expression of a single gene.

Then we will look at linear models for multiple genes using limma.


```r
library(ABC2018)
library(ggplot2)
theme_set(theme_minimal())
data("lmExamples")
names(lmExamples)
```

```
## [1] "twoGroup1"    "twoGroup2"    "interaction1" "interaction2"
## [5] "batchEffect"  "threeGroup"
```

# Linear models of a single gene

## A simple example

Let's consider the first and simplest dataset:


```r
lmExamples$twoGroup1
```

```
##    Group Expression
## 1   Ctrl   4.577016
## 2   Ctrl   3.450122
## 3   Ctrl   4.935571
## 4   Ctrl   5.270881
## 5   Ctrl   6.735284
## 6    Trt   6.735289
## 7    Trt   9.099471
## 8    Trt   7.863351
## 9    Trt   6.389413
## 10   Trt   7.637056
```

## A simple example

Let's plot the two groups:


```r
ggplot(lmExamples$twoGroup1, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

## A simple example

We might compare the means using a t-test:


```r
t.test(Expression~Group, data=lmExamples$twoGroup1, var.equal = TRUE)
```

```
## 
## 	Two Sample t-test
## 
## data:  Expression by Group
## t = -3.5746, df = 8, p-value = 0.007245
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -4.196891 -0.905391
## sample estimates:
## mean in group Ctrl  mean in group Trt 
##           4.993775           7.544916
```

## A simple example

In this simple example, lm() is identical to the t-test, since they are both based on minimizing the residual sum of squares (RSS):


```r
summary(lm(Expression~Group, data=lmExamples$twoGroup1))
```

```
## 
## Call:
## lm(formula = Expression ~ Group, data = lmExamples$twoGroup1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.54365 -0.71141  0.01697  0.30810  1.74151 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.9938     0.5046   9.896 9.18e-06 ***
## GroupTrt      2.5511     0.7137   3.575  0.00724 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.128 on 8 degrees of freedom
## Multiple R-squared:  0.615,	Adjusted R-squared:  0.5668 
## F-statistic: 12.78 on 1 and 8 DF,  p-value: 0.007245
```

What is the meaning of the coefficients?

## A simple example

`lm()` models the data as :

$$ Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i, i=1,\dots,n $$
This can be expressed using linear algebra terms as:

$$ \mathbf{Y}=\mathbf{X}\boldsymbol{\beta}+\boldsymbol{\varepsilon} $$

`lm()` calculate the coefficients in $\beta$ using Ordinary Least Squares (OLS):

First `lm()` build a design (or model) matrix, which is a presence-absence matrix of effects.


```r
mod <- model.matrix(~Group, data=lmExamples$twoGroup1)
```

## A simple example


```r
mod
```

```
##    (Intercept) GroupTrt
## 1            1        0
## 2            1        0
## 3            1        0
## 4            1        0
## 5            1        0
## 6            1        1
## 7            1        1
## 8            1        1
## 9            1        1
## 10           1        1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$Group
## [1] "contr.treatment"
```

## A simple example

Then `lm()` obtains OLS coefficients using: $$ \hat{\boldsymbol{\beta}} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y} $$


```r
solve(t(mod) %*% mod) %*% t(mod) %*% lmExamples$twoGroup1$Expression
```

```
##                 [,1]
## (Intercept) 4.993775
## GroupTrt    2.551141
```

See your favourite statistics text book for a proof of this!

Note: Actually `lm()` uses something called QR-decomposition for numeric stability, but the idea is the same.

## A simple example

`lm()` calculate p-values based on the assumption that noise is normally distributed.

In the previous example the difference was signficant. In the second examples it is not!

```r
summary(lm(Expression~Group, data=lmExamples$twoGroup2))
```

```
## 
## Call:
## lm(formula = Expression ~ Group, data = lmExamples$twoGroup2)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.9541 -0.7018 -0.2077  0.7866  1.0776 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7685     0.3720   12.82 1.29e-06 ***
## GroupTrt     -0.6417     0.5261   -1.22    0.257    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8318 on 8 degrees of freedom
## Multiple R-squared:  0.1568,	Adjusted R-squared:  0.05143 
## F-statistic: 1.488 on 1 and 8 DF,  p-value: 0.2573
```

## A simple example


```r
ggplot(lmExamples$twoGroup2, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

## Interaction design

The power of `lm()` over `t.test()` lies in the fact that `lm()` can estimate means in much more complicated models, where we have more than two means to estimate!

Consider for example a classic interaction experiment:


```r
lmExamples$interaction1
```

```
##    Group Time Expression
## 1   Ctrl   0h   6.263864
## 2   Ctrl   0h   5.250198
## 3   Ctrl   0h   5.258195
## 4   Ctrl   0h   6.785534
## 5   Ctrl   0h   3.780294
## 6    Trt   0h   6.759878
## 7    Trt   0h   5.941338
## 8    Trt   0h   7.419409
## 9    Trt   0h   6.729043
## 10   Trt   0h   6.368175
## 11  Ctrl   1h   7.771588
## 12  Ctrl   1h   9.178680
## 13  Ctrl   1h   7.733727
## 14  Ctrl   1h   8.528141
## 15  Ctrl   1h   6.231341
## 16   Trt   1h   9.510209
## 17   Trt   1h  11.410522
## 18   Trt   1h   8.924736
## 19   Trt   1h  10.292395
## 20   Trt   1h   9.793336
```

## Interaction design


```r
ggplot(lmExamples$interaction1, aes(x=paste0(Group, "-", Time), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

## Interaction design


```r
summary(lm(Expression~Group*Time, data=lmExamples$interaction1))
```

```
## 
## Call:
## lm(formula = Expression ~ Group * Time, data = lmExamples$interaction1)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.6873 -0.3256 -0.1360  0.6735  1.4243 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       5.4676     0.4312  12.681 9.21e-10 ***
## GroupTrt          1.1760     0.6097   1.929   0.0717 .  
## Time1h            2.4211     0.6097   3.971   0.0011 ** 
## GroupTrt:Time1h   0.9216     0.8623   1.069   0.3010    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9641 on 16 degrees of freedom
## Multiple R-squared:  0.7901,	Adjusted R-squared:  0.7508 
## F-statistic: 20.08 on 3 and 16 DF,  p-value: 1.134e-05
```

Is the interaction significant?

## Interaction design


```r
model.matrix(Expression~Group*Time, data=lmExamples$interaction1)
```

```
##    (Intercept) GroupTrt Time1h GroupTrt:Time1h
## 1            1        0      0               0
## 2            1        0      0               0
## 3            1        0      0               0
## 4            1        0      0               0
## 5            1        0      0               0
## 6            1        1      0               0
## 7            1        1      0               0
## 8            1        1      0               0
## 9            1        1      0               0
## 10           1        1      0               0
## 11           1        0      1               0
## 12           1        0      1               0
## 13           1        0      1               0
## 14           1        0      1               0
## 15           1        0      1               0
## 16           1        1      1               1
## 17           1        1      1               1
## 18           1        1      1               1
## 19           1        1      1               1
## 20           1        1      1               1
## attr(,"assign")
## [1] 0 1 2 3
## attr(,"contrasts")
## attr(,"contrasts")$Group
## [1] "contr.treatment"
## 
## attr(,"contrasts")$Time
## [1] "contr.treatment"
```

## Interaction design


```r
model.matrix(Expression~Group+Time+Group:Time, data=lmExamples$interaction1)
```

```
##    (Intercept) GroupTrt Time1h GroupTrt:Time1h
## 1            1        0      0               0
## 2            1        0      0               0
## 3            1        0      0               0
## 4            1        0      0               0
## 5            1        0      0               0
## 6            1        1      0               0
## 7            1        1      0               0
## 8            1        1      0               0
## 9            1        1      0               0
## 10           1        1      0               0
## 11           1        0      1               0
## 12           1        0      1               0
## 13           1        0      1               0
## 14           1        0      1               0
## 15           1        0      1               0
## 16           1        1      1               1
## 17           1        1      1               1
## 18           1        1      1               1
## 19           1        1      1               1
## 20           1        1      1               1
## attr(,"assign")
## [1] 0 1 2 3
## attr(,"contrasts")
## attr(,"contrasts")$Group
## [1] "contr.treatment"
## 
## attr(,"contrasts")$Time
## [1] "contr.treatment"
```

## Interaction design


```r
ggplot(lmExamples$interaction2, aes(x=paste0(Group, "-", Time), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

## Interaction design


```r
summary(lm(Expression~Group*Time, data=lmExamples$interaction2))
```

```
## 
## Call:
## lm(formula = Expression ~ Group * Time, data = lmExamples$interaction2)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.32644 -0.83516 -0.08706  0.37706  2.33671 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       4.6754     0.4793   9.754 3.88e-08 ***
## GroupTrt          2.4256     0.6779   3.578 0.002512 ** 
## Time1h            2.9263     0.6779   4.317 0.000532 ***
## GroupTrt:Time1h   4.3643     0.9587   4.553 0.000326 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.072 on 16 degrees of freedom
## Multiple R-squared:  0.9341,	Adjusted R-squared:  0.9217 
## F-statistic: 75.57 on 3 and 16 DF,  p-value: 1.155e-09
```

Is the interaction significant?

## Batch correction

Another common case is blocking or correcting for batches:


```r
lmExamples$batchEffect
```

```
##    Group Batch Expression
## 1   Ctrl     A   5.835183
## 2   Ctrl     A   5.589023
## 3   Ctrl     A   5.170094
## 4   Ctrl     A   4.562828
## 5   Ctrl     A   4.919792
## 6    Trt     A   6.072457
## 7    Trt     A   8.518867
## 8    Trt     A   7.085730
## 9    Trt     A   9.152795
## 10   Trt     A   5.634069
## 11  Ctrl     B  13.345906
## 12  Ctrl     B  13.565181
## 13  Ctrl     B  12.993902
## 14  Ctrl     B  12.561219
## 15  Ctrl     B  13.688891
## 16   Trt     B  14.500581
## 17   Trt     B  13.131385
## 18   Trt     B  15.743425
## 19   Trt     B  14.910825
## 20   Trt     B  16.081164
```

## Batch correction


```r
ggplot(lmExamples$batchEffect, aes(x=paste0(Group, "-", Batch), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

## Batch correction

Without modelling batch:

```r
summary(lm(Expression~Group, data=lmExamples$batchEffect))
```

```
## 
## Call:
## lm(formula = Expression ~ Group, data = lmExamples$batchEffect)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -5.449 -4.011  0.059  3.901  4.998 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    9.223      1.335   6.909 1.85e-06 ***
## GroupTrt       1.860      1.888   0.985    0.338    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.222 on 18 degrees of freedom
## Multiple R-squared:  0.05116,	Adjusted R-squared:  -0.001556 
## F-statistic: 0.9705 on 1 and 18 DF,  p-value: 0.3376
```

## Batch correction

With modelling batch:

```r
summary(lm(Expression~Group+Batch, data=lmExamples$batchEffect))
```

```
## 
## Call:
## lm(formula = Expression ~ Group + Batch, data = lmExamples$batchEffect)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.85083 -0.50149 -0.08485  0.52495  1.96875 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   5.3241     0.3845   13.85 1.09e-10 ***
## GroupTrt      1.8599     0.4439    4.19 0.000615 ***
## BatchB        7.7982     0.4439   17.57 2.47e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9927 on 17 degrees of freedom
## Multiple R-squared:  0.9505,	Adjusted R-squared:  0.9446 
## F-statistic: 163.1 on 2 and 17 DF,  p-value: 8.083e-12
```

## Batch correction


```r
model.matrix(~Group+Batch, data=lmExamples$batchEffect)
```

```
##    (Intercept) GroupTrt BatchB
## 1            1        0      0
## 2            1        0      0
## 3            1        0      0
## 4            1        0      0
## 5            1        0      0
## 6            1        1      0
## 7            1        1      0
## 8            1        1      0
## 9            1        1      0
## 10           1        1      0
## 11           1        0      1
## 12           1        0      1
## 13           1        0      1
## 14           1        0      1
## 15           1        0      1
## 16           1        1      1
## 17           1        1      1
## 18           1        1      1
## 19           1        1      1
## 20           1        1      1
## attr(,"assign")
## [1] 0 1 2
## attr(,"contrasts")
## attr(,"contrasts")$Group
## [1] "contr.treatment"
## 
## attr(,"contrasts")$Batch
## [1] "contr.treatment"
```

## Three Groups

As you can see, you can specify any complex design! For example a three group design:


```r
summary(lm(Expression~Group, data=lmExamples$threeGroup))
```

```
## 
## Call:
## lm(formula = Expression ~ Group, data = lmExamples$threeGroup)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.7032 -0.9613  0.2510  0.7385  2.5940 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   5.1988     0.6144   8.462 2.11e-06 ***
## GroupB        0.8047     0.8689   0.926 0.372622    
## GroupC        3.7632     0.8689   4.331 0.000977 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.374 on 12 degrees of freedom
## Multiple R-squared:  0.6342,	Adjusted R-squared:  0.5733 
## F-statistic:  10.4 on 2 and 12 DF,  p-value: 0.002395
```

How does the model matrix look?

## Three groups


```r
ggplot(lmExamples$threeGroup, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)

## Going further

The R formula interface is almost a programming language in itself!

Due to it's expressiveness, it's used by huge number of packages, making it well worth the effort to learn.

See this cheatsheet for an overview:

https://ww2.coastal.edu/kingw/statistics/R-tutorials/formulae.html

By default, the first alphabetical category is chosen as the Intercept. If you do not want this, you can use the `relevel` function to specify the reference level of the factor.

This "memory" of a factor to store it's reference level is why R insists on coercing characters to factors!

# Linear models of multiple genes

## Introduction

Now we have seen how linear models can be used a analyse the expression of a single gene.

In a genomic setting, we wish to analyze the expression of all genes. This means that we are applying the same model many times to each individual genes.

Given the low sample size of most genomics experiments, each linear model has relatively low power to detect differential expression.

We can get around this problem by assuming all models are somehow similar, and share information between them, i.e. share information between genes.

Here we show the simplest implementation of this using limma-trend.


```r
library(limma)
library(edgeR)
data("zebrafish")
```

## Normalizing the EM

Like previosly, we first normalize and log-transform the EM:


```r
# Trim
above_one <- rowSums(zebrafish$Expression > 1)
trimmed_em <- subset(zebrafish$Expression, above_one > 4)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(y=dge, log=TRUE)
```

## Setting up the model matrix.

We then set up of simple design matrix

```r
mod <- model.matrix(~gallein, data=zebrafish$Design)
mod
```

```
##       (Intercept) galleintreated
## Ctl1            1              0
## Ctl3            1              0
## Ctl5            1              0
## Trt9            1              1
## Trt11           1              1
## Trt13           1              1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$gallein
## [1] "contr.treatment"
```

## Gene-wise linear models

We then use limma to fit gene-wise linear models

```r
fit <- lmFit(EM, design=mod)
fit
```

```
## An object of class "MArrayLM"
## $coefficients
##                    (Intercept) galleintreated
## ENSDARG00000000001    2.272137    -0.47391987
## ENSDARG00000000002    3.456370    -0.14747423
## ENSDARG00000000018    2.466553     1.70405393
## ENSDARG00000000019    6.555697     1.17105985
## ENSDARG00000000068    1.156748     0.05061264
## 18101 more rows ...
## 
## $rank
## [1] 2
## 
## $assign
## [1] 0 1
## 
## $qr
## $qr
##       (Intercept) galleintreated
## Ctl1   -2.4494897     -1.2247449
## Ctl3    0.4082483      1.2247449
## Ctl5    0.4082483      0.2898979
## Trt9    0.4082483     -0.5265986
## Trt11   0.4082483     -0.5265986
## Trt13   0.4082483     -0.5265986
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$gallein
## [1] "contr.treatment"
## 
## 
## $qraux
## [1] 1.408248 1.289898
## 
## $pivot
## [1] 1 2
## 
## $tol
## [1] 1e-07
## 
## $rank
## [1] 2
## 
## 
## $df.residual
## [1] 4 4 4 4 4
## 18101 more elements ...
## 
## $sigma
## ENSDARG00000000001 ENSDARG00000000002 ENSDARG00000000018 
##          1.7402541          1.2458209          0.5093126 
## ENSDARG00000000019 ENSDARG00000000068 
##          0.9365106          1.2362995 
## 18101 more elements ...
## 
## $cov.coefficients
##                (Intercept) galleintreated
## (Intercept)      0.3333333     -0.3333333
## galleintreated  -0.3333333      0.6666667
## 
## $stdev.unscaled
##                    (Intercept) galleintreated
## ENSDARG00000000001   0.5773503      0.8164966
## ENSDARG00000000002   0.5773503      0.8164966
## ENSDARG00000000018   0.5773503      0.8164966
## ENSDARG00000000019   0.5773503      0.8164966
## ENSDARG00000000068   0.5773503      0.8164966
## 18101 more rows ...
## 
## $pivot
## [1] 1 2
## 
## $Amean
## ENSDARG00000000001 ENSDARG00000000002 ENSDARG00000000018 
##           2.035177           3.382633           3.318580 
## ENSDARG00000000019 ENSDARG00000000068 
##           7.141227           1.182054 
## 18101 more elements ...
## 
## $method
## [1] "ls"
## 
## $design
##       (Intercept) galleintreated
## Ctl1            1              0
## Ctl3            1              0
## Ctl5            1              0
## Trt9            1              1
## Trt11           1              1
## Trt13           1              1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$gallein
## [1] "contr.treatment"
```

## Empirical Bayes

Limma uses empirical bayes to _shrink_ t-statistics of genes toward the overall trend:

```r
eb <- eBayes(fit, trend=TRUE, robust=TRUE)
plotSA(eb)
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29-1.png)

## Empirical Bayes

We can see how the shrinage is working, by comparing the shrunken to the unshrunken t-statistics:


```r
tstat <- data.frame(raw=(fit$coef/fit$stdev.unscaled/fit$sigma)[,"galleintreated"],
										shrunken=eb$t[,"galleintreated"])
```

## Empirical Bayes


```r
ggplot(tstat, aes(x=raw, y=shrunken, color=raw-shrunken)) + 
	geom_point(alpha=0.33) +
	scale_color_distiller(palette = "Spectral") +
	geom_abline(linetype="dashed", alpha=0.75)
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31-1.png)

## Obtaining the shrunken estimates


```r
topTable(eb, coef="galleintreated")
```

```
##                        logFC    AveExpr         t      P.Value adj.P.Val
## ENSDARG00000089384  4.078938  6.2292934  5.066452 0.0004937548 0.9694505
## ENSDARG00000086686  4.951126  0.8939408  4.444266 0.0012486892 0.9694505
## ENSDARG00000088436  4.137011  1.4408195  4.274795 0.0016161282 0.9694505
## ENSDARG00000091744  3.914578  5.2309892  4.265627 0.0016640929 0.9694505
## ENSDARG00000002508 -6.893081 -1.2836386 -4.231247 0.0017567359 0.9694505
## ENSDARG00000079302 -5.808159  3.9069506 -4.141658 0.0020246322 0.9694505
## ENSDARG00000062055  2.674427  6.5573718  3.955166 0.0026842265 0.9694505
## ENSDARG00000028581  2.826668  2.4144751  3.886053 0.0030037756 0.9694505
## ENSDARG00000052073  2.714842  5.1794133  3.878474 0.0030525572 0.9694505
## ENSDARG00000090689  5.202294  1.9178292  3.862610 0.0031714945 0.9694505
##                            B
## ENSDARG00000089384 -4.557240
## ENSDARG00000086686 -4.560526
## ENSDARG00000088436 -4.561523
## ENSDARG00000091744 -4.561723
## ENSDARG00000002508 -4.561947
## ENSDARG00000079302 -4.562544
## ENSDARG00000062055 -4.563687
## ENSDARG00000028581 -4.564200
## ENSDARG00000052073 -4.564294
## ENSDARG00000090689 -4.564536
```

