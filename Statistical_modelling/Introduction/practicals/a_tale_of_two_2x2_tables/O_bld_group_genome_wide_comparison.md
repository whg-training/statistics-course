## Calibration / model checking for the O blood group finding

In class we analysed this table:

```
                      non-O blood group    O blood group
controls              3420                 3233
severe malaria cases  3925                 2738
```

which can be created in R as:
```
data = matrix(
    c(
        3420, 3233,
        3925, 2738
    ),
    ncol = 2,
    byrow = T,
    dimnames = list(
        c( "controls", "cases" ),
        c( "non-O", "O" )
    )
)
```

The odds ratio computed from this table is 0.74 - suggesting O blood group is protective.  And a test of this table using `chisq.test()` or `fisher.test()` shows that this is highly statistically significant.  That means:

* *If* the data were truly sampled from a binomial distribution with given frequency in each row
* and *if* the true odds ratio was zero
* *then* a table with such a large odds ratio would almost never occur.

These are approximations which don't really hold.  Samples are not really independent (e.g. population structure is itself a reflection of distant relationships); frequencies might vary across populations; in a finite population it's not really credible that the true effect size is really zero (though it might be very close).

We would therefore like to verify that there isn't something going on we haven't thought of.  I.e. is our null model really sensible?

One of the advantages of genome-wide analysis is that we have a lot of data to answer this with.
Suppose we compare our finding with all other 'similar' variants in the genome.  They have all been
genotyped on the same set of samples.  So they represent the same sampling biases (if any) and the same demography.  However, due to recombination they also represent many independent draws from the genealogical history of the sample.

The file `recessive_counts_matching_O_frequency.csv.gz` contains data for all the alleles in our
study that match the frequency of the O blood group mutation.

Load it:

```
X = read.csv(
    "recessive_counts_matching_O_frequency.csv.gz",
    header = T,      # Tell R our file has a row of column names
    as.is = T,       # Ask R to please not transmogrify any data
    check.names = F  # Also ask R not to mess around with column names 
)
```
Let's first sanity check I've computed frequencies right.  It's often best to first collect numerical data into a matrix - this is a container of numbers that is layed out nicely in memory for computation.  (Unlike a data frame which holds columns of data of different types).  The `as.matrix()` function does this:
```
data = as.matrix( X[, c( "controls:A", "controls:B", "cases:A", "cases:B" )])
# look at the data
head( data )

X$frequency = (data[,2] + data[,4]) / rowSums(data)

# plot it
hist(
    X$frequency,
    breaks = 5, # how many bars across the plot?
    xlab = "Recessive dosage frequency",
    xlim = c( 0.3, 0.7 ) 
)
abline( v = X$frequency[ which( X$rsid == 'rs8176719' )], col = 'red' )
```

Compute the log-odds ratio for every variant:
```
X$logOR = log( (data[,1] * data[,4]) / ( data[,2] * data[,3]) )
```

Let's look at the distribution of log odds ratios:
```
hist(
    X$logOR,
    breaks = 50,
    xlab = "log Odds ratio",
    main = "Alleles matching O blood group frequency",
    xlim = c( -0.5, 0.5 ),
    xaxt = 'n'  # Turn off the axis so we can draw our own
)
ticks = c( 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3 )
axis(
    1,       # on the bottom edge
    at = log(ticks),
    label = ticks
)

# Let's put an arrow where O blood group allele is:
arrows(
    x0 = X$logOR[ which( X$rsid == 'rs8176719' )],
    x1 = X$logOR[ which( X$rsid == 'rs8176719' )],
    y0 = 2000, y1 = 100,
    length = 0.1,
    col = 'red'
)
text(
    x = X$logOR[ which( X$rsid == 'rs8176719' )], y = 2000,
    pos = 3, # draw text above
    label = "O bld grp.",
    col = 'red'
)
```


Compute the standard error of the log OR (using the standard formula):
```
# for table
#   a b
#   c d
# the formula is sqrt( 1/a + 1/b + 1/c + 1/d )
X$SE = sqrt( rowSums( 1/data ))
```

Is the null model at all plausible here?  The null model says that all variation in the effect size (the log OR) is due to sampling variation.  It should be distributed as a N(0, SE^2).  (all the stnadard errors will be similar here because of the choice of variant).  To test how well the null model fits, we can plot standard-error-adjusted effect sizes ("Z-scores") against a gaussian distribution:
```
region = c( -5, 5 )
hist(
    X$logOR / X$SE,
    breaks = 100,
    xlab = "x score",
    xlim = region,
    ylim = c( 0, 0.5 ),
    main = "Alleles matching O blood group frequency",
    freq = FALSE # plot y axis values as a density, for comparison
)
x = seq( from = region[1], to = region[2], by = 0.01 )
points(
    x,
    dnorm( x, mean = 0, sd = 1 ),
    type = 'l',
    col = 'red'
)
```

It's not a bad fit!  There's a small gap at the top and slightly fatter tail to the real data.
A more sensitive plot is a quantile-quantile plot, where we plot ordered statistics against their expection.  We can do that here:
```
zscores = sort( X$logOR / X$SE )
expected = qnorm( (1:length(zscores))/(length(zscores)+1)) # Normal distribution quantiles
plot(
    expected,
    zscores,
    pch = 19 # nice round dots
)
abline( a = 0, b = 1, col = 'red' )
grid()
```
The distribution definitely has fatter tails than expected.

## Aside
(A: I've the qq plot in terms of z-scores, because I like thinking in 'effect size' space.  A more often seen qq-plot is based on a p-value as follows.  A traditional (Wald test) P-value compares
the z-score with the standard gaussian, as in our plot above:
```
X$pvalue = pnorm( -abs(X$logOR), sd = X$SE ) * 2
```
*Note*: the `abs()` and multiplying by two in the above are needed to capture both the left and right tail of the distribution.  We'd be interested in effects going in either direction.  (You could equally well use the squared z-score and `pchisq`.)  We can make a qq-plot by plotting p-values against ranked expected quantiles:
```
pvalues = -log10( sort( X$pvalue ) )
expected = -log10( (1:length(pvalues)  / (length(pvalues)+1)))
plot(
    expected,
    xlab = "Expected log10 P-value",
    pvalues,
    ylab = "Observed log10 P-value",
    pch = 19, # nice round dots
    xlim = c( 0, 20 ), # make it square
    ylim = c( 0, 20 )  # make it square
)
abline( a = 0, b = 1, lty = 2 )

# Let's highlight O blood group on the plkot:
w = which( pvalues == -log10( X$pvalue[ X$rsid == 'rs8176719' ]))
points(
    expected[w],
    pvalues[w],
    col = 'red'
)
text(
    expected[w],
    pvalues[w],
    pos = 2,             # put it to left
    adj = 1,             # right-justify
    label = "O bld grp"
)
```
Challenge: the kth point in this plot is -log10 of a kth order statistic.  Use this fact to add vertical error bars for each point?

## Summary

We have shown

1. That O blood has an estimated protective effect in these data.
2. That such a large estimated effect is unlikely *under the formal model assumptions* if the effect were really zero (hence the small P-value).
3. That other similar variants genome-wide do seem fairly compatible with the model of no effect.  This gives us confidence that our P-value (computed against the null) reflects something realistic. In particular, there isn't some large systematic problem with our use of the P-value that would invalidate our analysis.
4. However, genome-wide variants don't seem 100% compatible with the model assumptions.

Note that:

* point 1 is purely about our statistical model and the data we have observed
* point 2 is a comparison of our data with a *hypothetical infinite set of unobserved data* generated under the model assumptions.
* point 3 is a comparison of our data with other real data generated in the same samples.
* point 4 suggests there is still further calibration to perform.

Of course we also know that the O blood group mutation (rs8176719) is a deletion of protein-coding sequence that alters red cell surface antigens, making us a priori likely to think this might be involved (compared to, say, a randomly chosen genetic variant).

What could we do re: point 4?  Are there actually many variants across the genome having nontrivial true effect on malaria susceptibility?  Or is there some systematic confounding going on?  (This seems not impossible, due to the diverse set of samples data was drawn from.)

