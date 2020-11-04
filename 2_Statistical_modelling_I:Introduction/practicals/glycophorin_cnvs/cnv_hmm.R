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

# This function implements the log-sum-exp formula.
# If x = log(y) is a vector of values expressed in log space,
# this function computes log( sum(y) ) in a way that remains
# accurate even when sum(exp(x)) would over- or underflow.
log.sum.exp <- function( x, na.rm = FALSE ) {
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


# given:
# 1. a set of between-site distances
# 2. a lambda value (switch rate)
# 3. a prior probability of each copy number state
#
# compute.transitions() returns an array of transition log-probabilities
# This is arranged so that, if tp is the result then
# tp[n,i,j] is the log-probability of transitioning from state i to state j at position n
compute.transitions <- function(
    distances,      # distance between points
    lambda,         # 'switchiness'
    prior           # prior state log-probabilities 
) {
    # Given distance d, probability of switching
    # to a new state is 1 - e^-lambda d
    # (This starts at 0 when d = 0 and increases to 1 when d is large).
    # If switching we pick a new state from the prior.
    # (NB. this implies the chain is stationary.)
    # We work in log space so convert expressions accordingly
    D = length( prior )

    result = array(
        NA,
        dim = c(
            length( distances ),
            D,
            D
        )
    )
    for( n in 1:length(distances)) {
        distance = distances[n]
        scaled.prior = prior + log( 1 - exp( -lambda * distance ))
        result[n,,] = matrix(
            rep( scaled.prior, D ),
            byrow = T,
            nrow = D,
            ncol = D
        )
        # Probability of no switch is e^-lambda d
        diag( result[n,,] ) = rep( -lambda * distance, D )
    }
    # Note: (exponentiated) rows now sum to one because prior sums to 1.
    return( result )
}

# Given:
# 1: a set of coverage values (data)
# 2: modelled diploid mean and variance parameters for the sample
# 3: a set of copy numbers to consider
#
# compute.emissions() returns a matrix of emission probabilities.
#
# This is arranged so that if O is the return value then
# O[i,k] is the emission log-probability for the ith site and the kth copy number
compute.emissions <- function(
    data,               # vector of coverage values
    parameters = c(
        diploid.mean = 1,
        diploid.variance = 1
    ),
    copy.numbers       # vector of copy numbers
) {
    result = array(
        NA,
        dim = c(
            length(data),
            length(copy.numbers)
        ),
        dimnames = list(
            sprintf( "n=%d", 1:length(data )),
            sprintf( "cn=%d", copy.numbers )
        )
    )
    for( k in 1:length( copy.numbers )) {
        # multiply the mean and variance by the copy number.
        # However, there's always a bit of coverage due to spurious mapping
        # We deal with this by representing copy number 0 as 0.01.
        multiplier = max( copy.numbers[k], 0.01 ) / 2
        result[,k] = gaussian.ll(
            data,
            parameters = c(
                mean = as.numeric( parameters['diploid.mean'] * multiplier ),
                variance = as.numeric( parameters['diploid.variance'] * multiplier )
            )
        )
    }
    # If any data points were missing, we set the emission probability to 1
    # This means missing observations do not affect the Markov chain.
    result[is.na(result)] = 1
    return( result )
}


# forward.backward() implements (guess what?) the forward-backward algorithm
# based on arrays of emission and transition probabilities, and the prior
# First, the alpha and beta values (forward and backward algorithm probabilities) are computed.
# Then, the two are combined to compute state probabilities (gamma) at each site.
#
# These computations and notation are as in the Rabiner HMM tutorial:
# https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/tutorial%20on%20hmm%20and%20applications.pdf
#
# The only complication is that to avoid numerical over/underflow we work in log space,
# so have to change expressions accordingly.  In particular multiplications get converted
# to additions in log space.  These are best computed by the log-sum-exp formula.
#
forward.backward <- function(
    emissions,              # matrix of emission log-probabilities (LxD)
    transitions,
    prior                  # vector of prior state probabilities
) {
    # FORWARD ALGORITHM
    forward = function( emissions, transitions, prior ) {
        L = nrow(emissions)
        D = ncol(transitions)
        alpha = matrix( 0, nrow = L, ncol = D )
    
        # 1. Initialisation
        # alpha_1 gets initialised from prior and observations multiplied.
        # (Multiplied probabilities, but we are working in log space so add them)
        alpha[1,] = emissions[1,] + prior

        # alpha_n+1 is obtained by summing over states i at time n
        # In non-log-space this would be:
        #          ( sum_i alpha[n,i] * transitions[n,i,j] ) * emissions[n,j]
        # But in log space it is instead:
        #   alpha[n+1,j] = logsumexp_i(alpha[n,i] + transitions[n,i,j]) + emissions[n,j]
        for( n in 1:(L-1)) {
            for( j in 1:D ) {
                alpha[n+1,j] = log.sum.exp( alpha[n,] + transitions[n,,j] ) + emissions[n+1,j]
            }
        }
        return( alpha )
    }

    backward <- function( emissions, transitions ) {
        L = nrow(emissions)
        D = ncol(transitions)

        beta = matrix( 0, nrow = L, ncol = D )
    
        # Initialisation, initialise to prob = 1 (or log prob = 0)
        beta[L,] = 0
    
        # Usually
        # beta_n is obtained by summing over states j at time n+1
        #   sum_j transitions[n,i,j] * emissions[n+1,j] * beta[n+1,j]
        # but in log space it is instead:
        #   logsumexp_j( transitions[n,i,j] + emissions[n+1,j] + beta[n+1,j] )
    
        for( n in (L-1):1 ) {
            for( i in 1:D ) {
                beta[n,i] = log.sum.exp( transitions[n,i,] + emissions[n+1,] + beta[n+1,] )
            }
        }
        return( beta )
    }

    L = nrow(emissions)
    copy.numbers = 0:(length(prior)-1)

    # Posterior state probabilities
    alpha = forward( emissions, transitions, prior )
    beta = backward( emissions, transitions )

    # Posterior state probabilities
    gamma = alpha + beta - log.sum.exp( alpha[L,] )
    
    return(
        list(
            emissions = emissions,
            alpha = alpha,  # Forward log-probabilities
            beta = beta,    # Backward log-probabilities
            gamma = gamma # Marginal log-probabilities
        )
    ) ;

    # Sanity checks:
    # 1: rowSums( exp( gamma )) == 1, marginal probalities sum to one
    # 2: rowSums( exp( alpha + beta )) should all be equal (total observation probability)
}

# Implement one iteration of the CNV-finding HMM for a matrix of N samples and L sites (bins)
# Arguments:
# - `data` is a matrix of coverage values, with one site per row and one sample per column
# - `coverage.parameters` a dataframe containing modelled mean and variance for diploid copy number for each sample
# - `prior` contains state priors for the HMM
# - `lambda` is the HMM switch rate parameter, e,g, here lambda=1 means expect one switch per 1600bp bin
# - `site.multipliers` is a vector of values per site used to handle per-site variation in coverage
#
# Return value is a list with several members
# - `prior` is the copy number prior, as passed in
# - `marginal.log.probabilities` are the the log-posterior probabilities of each copy number given the HMM (i.e. gamma in the HMM)
# - `expected.copy.numbers` are the expected copy numbers at each site, given the posteriors.
# We also return the emission log-probabilities for downstream use.
cnv.hmm <- function( data, coverage.parameters, prior, lambda, site.multipliers ) {
    echo( "Running HMM for %d samples at %d sites...\n", ncol(data), nrow(data ))

    # prior is assumed to be for copy number states 0 ... K
    copy.numbers = 0:(length(prior)-1)
    log.prior = log( prior )

    #
    # This array reports the marginal probability of each
    # copy number state at each site for each sample
    result = list(
        prior = prior,
        marginal.log.probabilities = array(
            NA,
            dim = c( ncol(data), nrow(data), length( copy.numbers )),
            dimnames = list(
                colnames( data ),
                rownames( data ),
                sprintf( "cn=%d", copy.numbers )
            )
        ),
        expected.copy.numbers = array(
            NA,
            dim = c( nrow(data), ncol(data) ),
            dimnames = list(
                rownames( data ),
                colnames( data )
            )
        ),
        emission.log.probabilities = array(
            NA,
            dim = c( ncol(data), nrow(data), length( copy.numbers )),
            dimnames = list(
                colnames( data ),
                rownames( data ),
                sprintf( "cn=%d", copy.numbers )
            )
        ),
        total.log.probability = 0
    )
    transitions = compute.transitions( rep( 1, nrow(data)-1 ), lambda, log.prior )
    for( i in 1:ncol(data) ) {
        emissions = compute.emissions(
            data[,i] * site.multipliers,
            parameters = c(
                diploid.mean = coverage.parameters$means[i],
                diploid.variance = coverage.parameters$variances[i]
            ),
            copy.numbers
        )
        fb = forward.backward( emissions, transitions, log.prior )
        result$total.log.probability = result$total.log.probability + log.sum.exp( fb$alpha[1,] + fb$beta[1,] )
        result$marginal.log.probabilities[i,,] = fb$gamma
        result$expected.copy.numbers[,i] = exp(fb$gamma) %*% copy.numbers
        result$emission.log.probabilities[i,,] = emissions
        echo( "." )
        if( i %% 50 == 0 ) {
            echo( "\n" )
        }
    }
    echo( "\nok\n" )
    return( result ) ;
}

# Plot a matrix of copy numbers with a custom colour scheme
plot.copy.numbers <- function(
    copy.number,
    filename = NULL,
    title =  NULL,
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
sites = coverage[,1:3]
data = as.matrix( coverage[,4:ncol(coverage)] )
samples = colnames( data )

# We model coverage for diploid state as a gaussian with mean and variance
# We start by estimating these from the 1st 200 variants
coverage.parameters = data.frame(
    sample = colnames(data),
    means = sapply( 1:ncol( data ), function(i) { mean( data[1:200,i], na.rm = T ) }),
    variances = sapply( 1:ncol( data ), function(i) { var( data[1:200,i], na.rm = T ) })
)

# Run the HMM
result = cnv.hmm(
    # Binned coverage for each sample (columns) at each site (rows)
    data,
    # Mean and variance parameters for each sample
    coverage.parameters,
    # Prior on each copy number state
    # States are assumed to start at zero (no copies) and go up to K
    # Here we put 90% prior on diploid state
    prior = c( 0.02, 0.02, 0.9, 0.02, 0.02, 0.02 ),
    # lambda = switch rate per bin
    # Each bin is 1600bp, so 1/10 expects 1 switch every 16kb and so on
    lambda = 1/20,
    # site multipliers affect values across samples for each site
    # these are used to handle variation in coverage e.g. due to site-specific
    # mapping, sequence data rates, GC content etc.
    # Here we set these all to 1 for an initial run
    site.multipliers = rep( 1, nrow( data ))
)

plot.copy.numbers( result$expected.copy.number, "images/expected_posterior_copy_numbers_hmm.pdf", title = "Expected copy number (HMM model)" )
# Could now do MCMC over mean and variance parameters and over expected switch rate lambda

o = hclust( # hierarchical clustering
	    dist( t(result$expected.copy.number) ) # of Euclidean distance matrix between samples
)$order
plot.copy.numbers( result$expected.copy.number[,o], "images/expected_posterior_copy_numbers_hmm_clustered.pdf", title = "Expected copy number (HMM model, clustered)" )

