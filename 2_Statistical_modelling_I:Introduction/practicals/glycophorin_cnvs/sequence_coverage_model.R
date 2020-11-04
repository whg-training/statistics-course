echo <- function( message, ... ) {
    cat( sprintf( message, ... ))
}

gaussian.ll <- function(
    data,
    parameters = c(
        mean = 0,
        variance = 1
        )
) {
    return(
        dnorm( data, mean = parameters['mean'], sd = sqrt( parameters['variance'] ), log = TRUE )
    )
}

plot.copy.numbers <- function(
    copy.number,
    filename = NULL,
    title=  NULL,
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


# Load coverage data
# This data is sequence read coverage binned into 1600bp bins
# Bins are excluded (reported as NA) if they contain mostly 'unmappable' sites
# i.e. sites where reads can't be confidently aligned
coverage = read.delim( "glycophorin_binned_coverage.tsv.gz", header = T, as.is = T, sep = "\t" )
data = as.matrix( coverage[,4:ncol(coverage)] )
sites = coverage[,1:3]
samples = colnames( data )

# Here is what raw data looks like
image( data, x = 1:nrow(sites), y = 1:length(samples) )

# Let's model coverage for diploid state as gaussian with specified mean and variance
# We start by estimating these from the whole dataset
coverage.parameters = data.frame(
    sample = colnames(data),
    means = sapply( 1:ncol( data ), function(i) { mean( data[1:100,i], na.rm = T ) }),
    variances = sapply( 1:ncol( data ), function(i) { var( data[1:100,i], na.rm = T ) })
)


# If we think of the genome as 'time' and reads occur along the genome at a set rate
# then the number of reads starting at each position would be Poisson distributed
# with the given rate.
# Here we have binned reads into 1600bp bins, so scale up for approximation

pdf( file = "images/modelled_coverage.pdf", width = 9, height = 8 )
par( mar = c( 2, 2, 2, 1 ) ) # adjust margins to get more on plot
layout( matrix( 1:25, byrow = T, ncol = 5 ))
x = seq( from = 0, to = 20, by = 0.01 )
for( i in 1:25 ) {
	h = hist( data[,i], 25, freq = FALSE, main = colnames(data)[i] )
	points( x, dnorm( x, mean= coverage.parameters$means[i], sd = sqrt( coverage.parameters$variance[i] ) ), type = 'l', col = 'red' )
	grid()
}
dev.off()

# Let's make inference of copy number state just based on the (noisy) data at each site.
copy.numbers = 0:5 
copy.number.lls = array( NA, dim = c( nrow(data), ncol( data ), length( copy.numbers )), dimnames = list( sites$position, samples, sprintf( "cf=%d", copy.numbers )))
for( sample in 1:ncol(data)) {
	for( copy.number.index in 1:length(copy.numbers) ) {
		copy.number = copy.numbers[copy.number.index]
		mutiplier = max( copy.number, 0.01 ) # explained later
		copy.number.lls[,sample,copy.number.index] = gaussian.ll(
			data[,sample],
			parameters = c(
				mean = coverage.parameters$means[sample] * copy.number / 2,
				variance = coverage.parameters$variances[sample] * copy.number / 2
			)
		)
		# Handle missing data.
		# We will just treat missing data as having likelihood 1 (loglikelihood zero) here.
		copy.number.lls[ is.na(data[,sample]),sample, ] = 0
	}
}

# Find the maximum likelihood copy number state:
maximum.likelihood.state = array( NA, dim = c( nrow(data), ncol( data ) ), dimnames = list( sites$position, samples ))
for( variant in 1:nrow( data )) {
	for( sample in 1:ncol(data)) {
		w = which.max( copy.number.lls[variant,sample,] )
		stopifnot( !is.na(w) )
		maximum.likelihood.state[variant,sample] = copy.numbers[w]
	}
}

plot.copy.numbers( maximum.likelihood.state )

# That's daft, instead do a posterior probability
prior = c( `cn=0` = 0.02, `cn=1` = 0.02, `cn=2` = 0.9, `cn=3` = 0.02, `cn=4` = 0.02, `cn=5` = 0.02 ) # or choose your own
expected.posterior.state = array( NA, dim = c( nrow(data), ncol( data ) ), dimnames = list( sites$position, samples ))
for( variant in 1:nrow( data )) {
	for( sample in 1:ncol(data)) {
		log.unnormalised.posterior = copy.number.lls[variant,sample,] + log(prior)
		log.posterior = log.unnormalised.posterior - log.sum.exp( log.unnormalised.posterior )
		expected.posterior.state[variant,sample] = sum( exp(log.posterior) * copy.numbers )
		# or:
		# expected.posterior.state = sum( copy.number.lls[variant,sample,] + log(prior))
	}
}

plot.copy.numbers( expected.posterior.state )


# Look at them all together
pdf( "images/raw_binned_coverage.pdf", width = 6, height = 8 )
image( data, x = 1:nrow(sites), y = 1:length(samples), xlab = "Sites", ylab = "Samples", main = "Binned coverage" )
dev.off()

plot.copy.numbers( maximum.likelihood.state, "images/ml_copy_numbers.pdf", title = "maximum likelihood copy number (marginal)" )
plot.copy.numbers( expected.posterior.state, "images/expected_posterior_copy_numbers_nohmm.pdf", title = "Expected posterior copy number (marginal)" )

# And finally, cluster samples
o = hclust( # hierarchical clustering
	dist( t(expected.posterior.state) ) # of Euclidean distance matrix between samples
)$order
plot.copy.numbers( expected.posterior.state[,o], "images/expected_posterior_copy_numbers_nohmm_clustered.pdf", title = "Expected posterior copy number (marginal, clustered)" )
