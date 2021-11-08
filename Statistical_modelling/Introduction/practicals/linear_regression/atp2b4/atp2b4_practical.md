## Linear regression practical for *ATP2B4* expression data

In this part of the practical we will use linear regression to estimate the possible effect of a
genotype on gene expression values. The data is in the file `atp2b4_data.tsv` - load it now:

**Note.** I am using the [tidyverse](https://www.tidyverse.org) for these examples. If you don't
have it installed, either install it or use base R (e.g. `read.csv()`) instead.

```
library( tidyverse )
data = read_delim( 'atp2b4_data.tsv', delim = "\t" )
```

The data represents gene expression values (expressed as the mean per-site coverage of RNA-seq sequencing reads, divided by the total number of reads for the sample, in the column `normalised_covrage`) for the gene *ATP2B4*.  The data comes from [this paper](https://doi.org/10.1172/JCI94378) which reported that local genotypes - including at the SNP rs1419114 - are associated with expression changes.  Let's see if we can reconstruct this.

First let's rename the genotype column - too long!

```
colnames(data)[3] = "genotype"
```

**Question.** what does the data look like?  How many samples?  How many of each genotype?

### Fitting a regression model

Fit a linear model:
```
fit = lm( normalised_coverage ~ genotype, data = data )
coeffs = summary(fit)$coef
```

Let's plot:
```
plot(
  data$genotype,
  data$normalised_coverage,
  pch = 19,
  xlab = "rs1419114 genotype",
  ylab = "normalised_coverage",
  bty = 'n'
)
abline( coef = coeffs[,1], col = 'red', lwd = 2 )
abline( coef = coeffs[,1] - 1.96 * c( 0, coeffs[2,2] ), col = 'red', lty = 2, lwd = 1 )
abline( coef = coeffs[,1] + 1.96 * c( 0, coeffs[2,2] ), col = 'red', lty = 2, lwd = 1 )
```

And let's plot the residuals:
```
lm.qqplot( fit )
```
Hmm.  It's kind of ok.

If you look at the coefficients you'll see something like this:
```
> coeffs
             Estimate Std. Error  t value   Pr(>|t|)
(Intercept) 12.165760   4.785387 2.542273 0.01856574
genotype     4.395767   2.842943 1.546203 0.13632062
```

**Question.** What does that mean?

Hmm, this regression thinks there is a +ve effect of genotype on expression (4.4).  But the error interval is also wide (standard deviation of 2.8).

### getting the stage

However, this data was collected from two types of erythrocytes: fetal and adult.  Does adding that change the picture?
```
fit.with.stage = lm( normalised_coverage ~ stage + genotype, data = data )
coeffs.with.stage = summary(fit.with.stage)$coef

plot(
  data$genotype,
  data$normalised_coverage,
  pch = 19,
  xlab = "rs1419114 genotype",
  ylab = "normalised_coverage",
  bty = 'n'
)
abline( coef = coeffs.with.stage[c(1,3),1], col = 'red', lwd = 2 )

lm.qqplot( fit.with.stage )

# compare the coefficients:
print( coeffs.with.stage )
print( coeffs )
```

Suddenly things look different!

**Question.** Interpret these results.  What effect does `stage` have?  Which qq-plot looks better?

**Question.** which model (with stage or without) has a higher likelihood?  (Breaking my rule of not introducing P-values yet... you can do a 'likelihood ratio test':
```
# Compute likelihood ratio test statistic from log-likelihoods
log.lr = 2 * ( logLik( fit.with.stage ) - logLik( fit ))

# log.lr should be chi-square distributed if stage wasn't important
pchisq( log.lr, df = 1, lower.tail = F )
```
So stage substantially (and significantly) improves the model fit.

**Question.** If we put stage in, we get a strong genotype effect (3.5 standard deviations from zero).  If we don't put stage in, we get a much weaker effect. How do we know if we should include stage in the regression or not?

