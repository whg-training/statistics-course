## Model checking the malaria exposure finding

Consider the following table from <https://science.sciencemag.org/content/348/6235/711>:

```
                                        rs60822373 G allele        rs60822373 C allele
unexposed populations (Europeans)                      1965        1
malaria-exposed populations (Africans).                 707        17
```

It can be loaded in R like this:
```
theTable = matrix(
    c(
        1965, 1,
        707, 17
    ),
    ncol = 2,
    byrow = T,
    dimnames = list(
        c( "non-exposed", "exposed" ),
        c( "G", "C" )
    )
)
```

The odds ratio computed from this table is 47 - implying the frequencies of this variant are very different in European and African populations.
The Egan et al paper suggests this is evidence for malaria-driven natural selection of this allele.




How could we test that?  Well one way is 


And a test of this table ( `chisq.test()` or `fisher.test()` shows that this is highly statistically significant.  That means:

* *If* the data were truly sampled from a binomial distribution with given frequency in each row.
* and *if* the true odds ratio was zero.
* *then* a table with such a large odds ratio would almost never occur.

(*NB.* The assumptions of Fisher's exact test and of the chi-squared test, which corresponds to 'binomial sampling in rows', are slightly different.
Fisher's exact test is actually conditional on both row and column sums, which reduces things to one parameter.  This doesn't matter much for most tables in practice because the row and column sums are not usually very informative about the odds ratio itself.)

Does this table provide evidence that the malaria-exposed population allele frequency is higher because of natural selection due to malaria?  (That's the implicit message of this part of the original paper: [https://science.sciencemag.org/content/348/6235/711](Egan et al Science 2015) ).

One of the advantages of genome-wide analysis is that we have a lot of data to work with - and this can be used for model checking. Suppose
we compare our finding with all other 'similar' variants in the genome. They have all been
genotyped on the same set of samples, and as such they represent the same sampling biases (if any) and the same demography. However, due to recombination they also represent many independent draws from the genealogical history of the sample.

The file `allele_counts_with_exposed_count\=17.csv.gz` contains data for all the alleles in the 1000 Genomes populations studied, that have the same allele count of 17 in the exposed populations. Load it:

```
X = read.csv(
    "allele_counts_with_exposed_count=17.csv.gz",
    header = T,      # Tell R our file has a row of column names
    as.is = T,       # Ask R to please not transmogrify any data
    check.names = F  # Also ask R not to mess around with column names 
)

# Take a look
head(X)
```
Let's first sanity check I've computed the data right.
```
range( X$`exposed:B` )
```

It's often best to first collect numerical data into a matrix - this is a container of numbers that is layed out nicely in memory for computation.  (Unlike a data frame which holds columns of data of different types).  The `as.matrix()` function does this:
```
data = as.matrix( X[, c( "unexposed:A", "unexposed:B", "exposed:A", "exposed:B" )])
# Take a look
head( data )
```

Let's compute the frequency of each variant:
```
X$frequency = (data[,2] + data[,4]) / rowSums(data)

# plot it
hist(
    X$frequency,
    breaks = 50, # how many bars across the plot?
    xlab = "Allele frequency"
)
abline( v = X$frequency[ which( X$rsid == 'rs60822373' )], col = 'red' )
```

In this dataset we have conditioned on the exposed population frequency, so it's proabably more sensible to look at the frequency in non-exposed populations:

```
hist(
    X$`unexposed:B`,
    breaks = 50,
    xlab = "log Odds ratio",
    main = "Alleles with count = 17 in malaria-exposed populations"
)

# The value for our SNP is 1, so let's replot clamped to 50 to focus better on that area
hist(
    pmin( X$`unexposed:B`, 50 ),
    breaks = 50,
    xlab = "Count of allele in non-exposed pops",
    main = "Alleles with count = 17 in malaria-exposed populations",
    #xlim = c( -0.5, 0.5 ),
    #xaxt = 'n'  # Turn off the axis so we can draw our own
)
# Let's put a line where rs60822373 is:
abline(
    v = X$`unexposed:B`[ which( X$rsid == 'rs60822373' )],
    col = 'red'
)
text(
    x = X$`unexposed:B`[ which( X$rsid == 'rs60822373' )],
    y = 50000,
    adj = 0,
    pos = 4, # draw text to right
    label = "rs60822373",
    col = 'red'
)
```

Sure enough there are a bunch of alleles at higher counts in non-exposed populations.  But there is also a big spike at zero. 

We could also mesaure differences directly using the log-odds ratio:
```
X$logOR = log( ((data[,1]) * (data[,4])) / ( (data[,2]) * (data[,3])) )
```
Unfortunately if you do this you'll find many are infinite (due to a 0 in the table):
```
length( which( X$logOR == Inf ))
```

We will therefore compute the log OR by adding a dummy count of  1 to each cell.  This amounts to adding prior information.  This amounts to adding a prior.  (If you want to see what that prior looks like, plot the `table.ll.reparameterised()` function on a table full of 1's.)
```
X$logOR = log( ((data[,1]+1) * (data[,4]+1)) / ( (data[,2]+1) * (data[,3]+1)) )
```

Let's look at the distribution of the logOR:

```
# Compute log OR with additional prior count of 1 in each cell to avoid infinities.

hist(
    X$logOR,
    breaks = 50,
    xlab = "log Odds ratio",
    main = "Alleles with count = 17 in malaria-exposed populations"
)
abline(
    v = X$logOR[ which( X$rsid == 'rs60822373' )],
    col = 'red'
)
```
(Because the 3rd and 4th cells of our table are constant here, this has exactly the same information as above but expressed as a log OR).

We originally computed a P-value which indicated that the table was very unlikely to arise by chance under the null.  But now we find that many other variants have at least as extreme counts! We can use our data to get an empirical P-value:
```
empirical.P = length( which( X$`unexposed:B` <= 1 )) / nrow(X)
```

The answer is `0.49`.  So about half the SNPs in the genome with this low count in malaria-exposed populations, have equally low frequencies in non-exposed populations.  They can't all be under selection can they?

The null model says that all variation in the effect size (the log OR) is due to sampling variation true log OR = 0).  It should be distributed as a N(0, SE^2).  (all the standard errors will be similar here because of the choice of variant).  To test how well the null model fits, we can plot standard-error-adjusted effect sizes ("Z-scores") against a gaussian distribution.  First compute the standard error:

```
# for table
#   a b
#   c d
# the formula is sqrt( 1/a + 1/b + 1/c + 1/d )
X$SE = sqrt( rowSums( 1/(data+1) ))
```
*Warning*: this standard error is not quite correct due to our use of a prior additional count of 1 to each cell, which is fixed and not sampled.  This won't matter for this exercise.

```
region = c( -10, 10 )
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
The null model, at least as expressed through the log OR and its standard error, is not even close to being true.

For completeness / comparison, a more sensitive plot that we don't really need here, is a quantile-quantile plot, which plots ordered statistics against their expection.  We can do that here:
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
If the expected and true distributions were similar the points would lie near the line.

## Interpretation

Even though this table and the one for O blood group looked superficially similar, they differ in a couple of ways.  

One is that the small counts in some cells of this table mean we are at the boundary of where traditional statistics based on asymptotic assumptions holds.  This was the problem with the infinite log ORs and means we should be careful about just rolling out standard tools.  (The original authors used 'Fishers Exact Test' which does deal with this issue, although it introduces additional assumptions.)

Another more serious issue is that the assumed null model - that frequencies should be similar between the two groups - is not even close to being correct for most variants in the genome.  This is a serious issue which should cause us to reevaluate our assumptions.

