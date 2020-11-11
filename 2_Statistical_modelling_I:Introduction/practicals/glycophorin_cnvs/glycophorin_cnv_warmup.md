# Hidden Markov Model example warmup #
Welcome to the warm-up for the Hidden Markov Models session this week.

Data for this practical can be found in the `glycophorin_binned_coverage.tsv.gz` file.
And this practical produces plots - I've precomputed some of these and stored in the `images/` folder so you can see what they should look like.

The practical studies a region of the human genome on chromosome 4 that is known to contain large *copy number variants*.
These are genetic variants that occur (often) due to mistakes in DNA replication.  They lead to deletions, duplications of large tracts of DNA.
In general CNVs are extremely important determinants of human disease - the ones in this practical delete or duplicate whole genes!
But they are often also hard to access because they lie in 'difficult' regions of the genome, where paralogy (caused by ancestral duplications of DNA) makes analysis difficult.

In this practical the plan is to try to genotype some CNVs in one such region - the region containing GYPA, which
encodes one of the most abundant red cell surface proteins.  

To do this we will look at *sequence coverage data* (i.e. how many reads aligned to each genomic location from
short-read sequence data) and look for variation in coverage that might indicate loss or gain of DNA copies.
To make this simple, we have grouped the genome into consecutve 1600bp bins and we work with *mean coverage in each bin*.

## First look at the data

First of all let's load and look at the data - we'll use R for this practical (but if you are expert you are welcome to
explore other methods).

```R
data = read.delim( "glycophorin_binned_coverage.tsv.gz", header = T, as.is = T, sep = "\t" )
View(data)
```

The data has position information in the first few columns, and samples in columns 4 onwards.
Let's split these things out for easier handling.  We'll call the actual data `X`:
```R
X = as.matrix( data[,4:ncol(data)] )
sites = data[,1:3]
samples = colnames( X )
```

How could we plot this data?  One way is just to make a heatmap:

```
image( X, x = 1:nrow(sites), y = 1:length(samples) )
```

This is not all that edifying but a few features are evident. Clearly samples vary in coverage (some rows are more red
than othgers). Also, some sites seem to have more coverage than others (some columns have more coverage than others).
Also some sites have missing data! The reason for this is that the region contains paralogous gene copies which make
read mapping difficult - we have excluded bins where mappability was poor.

If you stare hard between bins 300-400, you'll see some samples seem to have something going on (long bands of more
yellow or more red).  This is actually the signal we want to extract, and to do that we need a model.

## Using an empirical model to handle variation in the data

Let's look at how coverage actually looks across sites, for the first few samples:
```R
par( mar = c( 2, 2, 2, 1 ) )                                    # this line adjust margins to get more on plot
layout( matrix( 1:25, byrow = T, ncol = 5 ))                    # put multiple panes on one plot
for( i in 1:25 ) {
    h = hist(
        X[,i],                      # coverage values for individual i
        breaks = 25,                # Number of bars to plot
        freq = FALSE,               # This scales the y axis so density sums to 1, i.e. empirical distribution
        main = samples[i]           # This is the plot title
    )
    grid()
}
```

The data for each sample looks kind of uni-modal and sort of symmetrical-ish in general, with a few bumps. Of course
any CNVs might affect that, and that could explain some of the bumps - that's what we want to find out.

For a practical approach we will assume that this binned coverage follows a Gaussian distribution. And for a first
guess, we will fit the gaussian using all the data in the 1st 100 bins (which, as we will see, happen not to contain
many CNVs). So let's go ahead and compute the parameters of these gaussians now:

```
coverage.parameters = data.frame(
    sample = samples,
    mean = sapply( 1:ncol( X ), function(i) { mean( X[1:100,i], na.rm = T ) }),
    variance = sapply( 1:ncol( X ), function(i) { var( X[1:100,i], na.rm = T ) })
)
```

We're going to use a Gaussian likelihood function to model sequence coverage data. Here's the function we wrote in the
Statistical modelling session (here implemented using `dnorm`).
```R
gaussian.ll <- function(
    data,
    params = list(
        mean = 0,
        variance = 1
    )
) {
    return(
        dnorm( data, mean = params$mean, sd = sqrt( params$variance ), log = TRUE )
    )
}
```

Let's plot data again and see how the fit looks:
```
par( mar = c( 2, 2, 2, 1 ) )                                    # this line adjust margins to get more on plot
layout( matrix( 1:25, byrow = T, ncol = 5 ))
x = seq( from = 0, to = 20, by = 0.01 )
for( i in 1:25 ) {
    h = hist(
        X[,i],                      # coverage values for individual i
        breaks = 25,                # Number of bars to plot
        freq = FALSE,               # This scales the y axis so density sums to 1, i.e. empirical distribution
        main = samples[i]           # This is the plot title
    )
    points(
        x,
        exp(gaussian.ll(
            data = x,
            params = list(
                mean = coverage.parameters$mean[i],
                variance = coverage.parameters$variance[i]
            )
        )),
        type = 'l',
        col = 'red'
    )
    grid()
}
```

The fit looks ... ok.  Not perfect, but not bad, depending on the sample.

*Exercise* feel free to look at other samples - you could also check model fit using a qq plot.

Clearly this model is not perfect but it might be enough to work with - for the purposes of this practical we'll go with it.

## Modelling copy number variation

Our general model of the effect of CNVs on coverage is that coverage of a site with copy number c should be (c/2) times
as large as the copy number at a diploid site - plus noise. If we think of sequence reads as being generated from each
copy independently, a bit of thought shows that the same relationship happens for the variance as well. So our model
for copy number c would be the *transformed Gaussian N( c * mean, c * variance )*.

This is good because it provides the key link from copy number to coverage, with noise, that might allow us to infer
copy numbers.

Let's try that now.  For simplicity:

- We will assume only copy numbers of 0 to 5 are possible (i.e. from a homozygous deletion up to a possible three extra copies).  Increase this if you want to!
- We'll put the results in a giant multidimensional array of sites x samples x copy number states.

We start by computing the per-site likelihood of each the data given each copy number.
Importantly as in the stats modelling session *we will work throughout in log space* to avoid any numerical issues.
This makes some computation more difficult, but 

```
copy.numbers = 0:5 
copy.number.lls = array(
    NA,                                 # fill with missing data to start
    dim = c(                            # dimensions
        nrow(sites),
        length(samples),
        length( copy.numbers )
    ),
    dimnames = list(                    # A good feature of R is you can name everything
        sites$position,
        samples,
        sprintf( "cf=%d", copy.numbers )
    )
)
```

Let's compute these lls:
```
for( sample in 1:length(samples) ) {
    for( i in 1:length(copy.numbers) ) {
        copy.number = copy.numbers[i]
        mutiplier = max( copy.number, 0.01 ) # explained later
        copy.number.lls[,sample,i] = gaussian.ll(
            X[,sample],
            params = list(
                mean = coverage.parameters$mean[sample] * copy.number / 2,
                variance = coverage.parameters$variance[sample] * copy.number / 2
            )
        )
        # Handle the missing data sites.
        # We will just treat missing data as having likelihood 1 (loglikelihood zero) here.
        copy.number.lls[ is.na(X[,sample]),sample, ] = 0
    }
}
```

Look at the top left of the output:
```
print( copy.number.lls[1:10,1:10,]
```

With any luck `copy.number.lls` captures some of the important signal in the data.  But how to plot it given it's 3-dimensional?  Here are two ways we could do it:

1. for each sample and site, we could take the maximum likelihood copy number call and plot that.
2. That's daft, as we're pretty sure most people are diploid at most sites. Instead let's apply our bayesian reasoning, and compute a posterior estimate of the copy number state, under a prior distribution that puts most weight on the diploid state.

The rest of this practical does both these things.  We will use this bit of code to actually make the plots:

```R
plot.copy.numbers <- function(
    copy.number,
    filename = NULL,
    title = NULL,
    palette = c( "darkorange2", "darkgoldenrod1", "forestgreen", "deepskyblue2", "deepskyblue4", "blueviolet" )
) { 
    if( !is.null( filename )) {
        pdf( file = filename, width = 6, height = 8 )
    } 
    par( mar = c( 4, 4, 2 + 2 *!is.null(title), 2 ) + 0.1 )
    image( copy.number, x = 1:nrow(copy.number), y = 1:ncol(copy.number), xlab = "sites", ylab = "samples", col = palette, zlim = c( 0, length(palette)-1 ), main = title )
    legend( "topleft", pch = 22, col = 'black', pt.bg = palette, legend = c( "hom deletion", "het deletion", "normal", "1 extra copy", "2 extra copies", "3 extra copies" ), bg = "white" )
    if( !is.null( filename )) {
        dev.off()
    }
}
```

## Copy number inference using MLEs

Ok - let's compute the maximum likelihood copy number state for each individual at each site:
```
maximum.likelihood.state = array( NA, dim = c( nrow(X), ncol( X ) ), dimnames = list( sites$position, samples ))
for( variant in 1:nrow( X )) {
    for( sample in 1:ncol(X)) {
        w = which.max( copy.number.lls[variant,sample,] )
        stopifnot( !is.na(w) )
        maximum.likelihood.state[variant,sample] = copy.numbers[w]
    }
}

plot.copy.numbers( maximum.likelihood.state, title = "Maximum likelihood copy number" )
```

That plot definitely looks cleaner!  Still lots of noise though, can we do better?

## Bayesian copy number inference

Taking the maximum likelihood estimate is daft though.  We have prior information: most people will be diploid at most sites with a few CNVs sprinkled in.  Let's build that knowledge in by applying our bayesian reasoning.  Still working seperately for each sample and each site, we will take a prior distribution that puts most weight on diploid state, and estiamte copy number state by taking posterior expectations.

```R
# our prior - note this should sum to one
prior = c(
    `cn=0` = 0.02,
    `cn=1` = 0.02,
    `cn=2` = 0.9,       # I put 90% weight on diploid state - you can try different values
    `cn=3` = 0.02,
    `cn=4` = 0.02, 
    `cn=5` = 0.02
) 
```

Here is another big array to put our results in:
```R
expected.posterior.state = array(
    NA,
    dim = c( nrow(sites), length(samples) ),
    dimnames = list(
        sites$position,
        samples
    )
)
```

Now let's compute it.

But wait! Here a technical problem occurs.  To compute the normalisation factor in Bayes rule we need to sum
over the possible copy numbers.  (Just like with the fair/unfair dice example fom the stats modelling session.)
But we are working in log space, so all values are log probabilities.
One way to do this would be to write `sum(exp(log.probs))` (and then take logs again to get back to log space).  However, a more numerically stable computation is to use the "log-sum-exp" formula.
This computes the same thing but avoids numerical errors by processing the largest value seperately.
Here is a fairly robust implementation:
```R
log.sum.exp <- function(
    x,
    na.rm = FALSE               # This mimics the na.rm argument of base R's sum() function
) {
    result = NA
    z = max( x )
    if( !is.na( z )) {
        terms = x-z
        if( na.rm ) {
            w = which( !is.na( terms ) )
        } else {
            w = 1:length( terms )
        }
        if( length( w ) > 0 ) {
            result = z + log( sum( exp( terms[w] ) ))
        }
    }
    return( result )
}
```

Now let's compute the posteriors:
```R
for( variant in 1:nrow( sites )) {
    for( sample in 1:length(samples)) {

        # These three lines implement Bayes rule, but working in log space
        log.unnormalised.posterior = copy.number.lls[variant,sample,] + log(prior)
        log.normalisation.constant = log.sum.exp( log.unnormalised.posterior )
        log.posterior = log.unnormalised.posterior - log.normalisation.constant
        
        # Now compute expected copy number given the posterior distribution
        expected.posterior.state[variant,sample] = sum( exp(log.posterior) * copy.numbers )
    }
}

plot.copy.numbers( expected.posterior.state, title = "Expected posterior copy number" )
```

Compare this with the maximum likelihood copy number calls above - the posterior version is much cleaner.

It's clear now there is something goin on in that region - a number of samples have clear runs of deleted or sometimes apparently duplicated bins.

To draw out the signal we can cluster samples - here using hierarchical clustering:

```
clustered.order = hclust(
    # we will cluster based on a distances between expected posterior state vectors
    dist( t(expected.posterior.state) )             
)$order

plot.copy.numbers( expected.posterior.state[,clustered.order], title = "Expected posterior copy number (clustered)" )
```

