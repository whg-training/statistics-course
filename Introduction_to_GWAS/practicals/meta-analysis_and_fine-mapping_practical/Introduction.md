[Up to the table of contents](README.md) - [Forward to the meta-analysis section](Meta-analysis.md)

## Meta-analysis and fine-mapping practical

So you've run your GWAS and you've found some signals.  Now what?

In this practical we will implement two steps that both try to help narrow down association signals
to the possible causal variants. In the first of these, we will take data from two studies and
*meta-analyse* it to combined the evidence. Then, we will use a fine-mapping method to try to identify a variant or variants that have the most evidence.

**Note**. Why can't we just use the P-values from the study directly, and pick the smallest?  Well maybe that would work, but there are two reasons to try something else:

* There might be more than one causal variant.
* The P-values reflect a mixture of *effect size* and *allele frequency*.  They might not be the best way to identify the underlying variants.

#### Getting the data

The data for this practical is in two files, `study1.z` and `study2.z`.  Load them now and take a look:

**Note.** I am using the [tidyverse](https://www.tidyverse.org) for these examples.  If you don't have that or don't want to get it, you can use base R functions like `read.table()` instead.

```
library(tidyverse)
study1 = read_delim( "study1.z", delim = " " )
study2 = read_delim( "study2.z", delim = " " )

```

Look at the data using `head()` or `View()`.

**Question.** What's in the data?  How many SNPs?

**Note.** For the analysis to work, the files will need to have the same SNPs in the same order.  Check this now!  (For example, you could check that `length( which( study1$rsid != study2$rsid ) )` is zero.)

The data for each study consists of:

* the rsid, chromosome, position and alleles of each SNP
* the minor allele frequency (`maf`) of each SNP
* and the GWAS / regression summary statistics, i.e. the maximum likelihood estimate *β* and its standard error `s`, for the association test of the SNP against the phenotype.

The data is for a region of the genome (overing *FUT2*) where an association has been found.

The data doesn't contain any P-values or Bayes factors. However, if you joined the regression
practical earlier you will know how to compute them:

```
  study1$P = pnorm( -abs(study1$beta), sd = study1$se ) * 2
```

and similarly for study2. The expression above computes the mass under the two tails with effect >=
*β*, of the gaussian distribution with standard deviation = the standard error. I.e. it is a
P-value using the effect size estimate as a test statistic.

**Note.** We could also compute a Bayes factor.  For example, suppose we believe in relatively small effects, so choose a *N(0,0.2<sup>2</sup>)* prior.  The calculation is:

```
study1$log10_BF = log10(
  pnorm( study1$beta, mean = 0, sd = sqrt( study1$se^2 + 0.2 ) ) /
  pnorm( study1$beta, mean = 0, sd = study1$se )
)
```

See [more here]().

**Question.** How much evidence is there in the two studies? Are any of the most-associated SNPs
shared between the two studies?

When you're ready, [start meta-analysing](Meta-analysis.md).



