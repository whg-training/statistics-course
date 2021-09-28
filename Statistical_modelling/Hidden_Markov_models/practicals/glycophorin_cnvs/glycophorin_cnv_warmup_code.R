# This file contains just the needed code portions 
# from glycophorin_cnv_warmup.md

# Log-sum-exp formula
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

# Gaussian ll
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

# Function to compute copy numbers
plot.copy.numbers <- function(
    copy.number,
    filename = NULL,
    title = NULL,
    palette = c( "darkorange2", "darkgoldenrod1", "forestgreen", "deepskyblue2", "deepskyblue4", "blueviolet" )
) { 
    if( !is.null( filename )) {
        pdf( file = filename, width = 6, height = 8 )
    } 
    par(mfrow=c(1,1))
    par( mar = c( 4, 4, 2 + 2 *!is.null(title), 2 ) + 0.1 )
    image( copy.number, x = 1:nrow(copy.number), y = 1:ncol(copy.number), xlab = "sites", ylab = "samples", col = palette, zlim = c( 0, length(palette)-1 ), main = title )
    legend( "topleft", pch = 22, col = 'black', pt.bg = palette, legend = c( "hom deletion", "het deletion", "normal", "1 extra copy", "2 extra copies", "3 extra copies" ), bg = "white" )
    if( !is.null( filename )) {
        dev.off()
    }
}

# Load data
data = read.delim( "glycophorin_binned_coverage.tsv.gz", header = T, as.is = T, sep = "\t" )
View(data)

# Split data
sites = data[,1:3]
X = as.matrix( data[,4:ncol(data)] )
rownames(X) = sites$position
samples = colnames( X )

# Compute coverage parameters
coverage.parameters = data.frame(
    sample = samples,
    mean = sapply( 1:ncol( X ), function(i) { mean( X[1:100,i], na.rm = T ) }),
    variance = sapply( 1:ncol( X ), function(i) { var( X[1:100,i], na.rm = T ) })
)

# Compute copy number lls
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

for( sample in 1:length(samples) ) {
    for( i in 1:length(copy.numbers) ) {
        copy.number = copy.numbers[i]
        mean.multiplier = copy.number/2
        variance.multiplier = max( copy.number/2, 0.01 ) # this to be error tolerant, as explained below
        copy.number.lls[,sample,i] = gaussian.ll(
            X[,sample],
            params = list(
                mean = coverage.parameters$mean[sample] * mean.multiplier,
                variance = coverage.parameters$variance[sample] * variance.multiplier
            )
        )
        # Handle the missing data sites.
        # We will just treat missing data as having likelihood 1 (loglikelihood zero) here.
        copy.number.lls[ is.na(X[,sample]),sample, ] = 0
    }
}

# Compute maximum likelihood state
maximum.likelihood.state = array( NA, dim = c( nrow(X), ncol( X ) ), dimnames = list( sites$position, samples ))
for( variant in 1:nrow( X )) {
    for( sample in 1:ncol(X)) {
        w = which.max( copy.number.lls[variant,sample,] )
        stopifnot( !is.na(w) )
        maximum.likelihood.state[variant,sample] = copy.numbers[w]
    }
}

# Specify our prior - note this should sum to one
prior = c(
    `cn=0` = 0.02,
    `cn=1` = 0.02,
    `cn=2` = 0.9,       # I put 90% weight on diploid state...
    `cn=3` = 0.02,      # ...and 2% weight on everything else - feel free to try different values (make it sum to 1)
    `cn=4` = 0.02, 
    `cn=5` = 0.02
) 
```

# Compute expected posterior state
expected.posterior.state = array(
    NA,
    dim = c( nrow(sites), length(samples) ),
    dimnames = list(
        sites$position,
        samples
    )
)

for( variant in 1:nrow( sites )) {
    for( sample in 1:length(samples)) {

        # These three lines implement Bayes rule, but working in log space
        log.unnormalised.posterior = copy.number.lls[variant,sample,] + log(prior)
        log.normalisation.constant = log.sum.exp( log.unnormalised.posterior )      # equivalent to: log( sum( exp( log.unnormalised.posterior )))
        log.posterior = log.unnormalised.posterior - log.normalisation.constant
        
        # Now compute expected copy number given the posterior distribution
        expected.posterior.state[variant,sample] = sum( exp(log.posterior) * copy.numbers )
    }
}

plot.copy.numbers( expected.posterior.state, title = "Expected posterior copy number" )
