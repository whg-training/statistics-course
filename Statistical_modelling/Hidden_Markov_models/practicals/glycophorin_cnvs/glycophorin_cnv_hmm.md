# Hidden Markov Model CNV calling example - continued #

** !! Note !! ** To run this code, you first need to have computed a set of log-likelihoods for each sample at each site in the relevant dataset.  The code to do this was explained in the [warmup tutorial](glycophorin_cnv_warmup.md), which you can also find in this directory.

If you haven't run it in this R session, please run this now:
```R
    source( 'glycophorin_cnv_warmup_code.R' )
```

You should now have the data (`X`, `sites`, `samples`) loaded, and you should have a 532 x 200 x 6 array called `copy.number.lls` loaded.  Check by writing:
```R
    View( copy.number.lls )
```

You should see a bunch of (mostly large and negative) numbers.

## The forward-backward algorithm

The forward.backward() function implements (guess what?) the HMM forward-backward algorithm based on arrays of emission
and transition probabilities, and the prior. First, the alpha and beta values (forward and backward algorithm
probabilities) are computed. Then, the two are combined to compute state probabilities (gamma) at each site.

These computations and notation are all as in the Rabiner HMM tutorial:
<https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/tutorial%20on%20hmm%20and%20applications.pdf>

The only complication is that to avoid numerical over/underflow we work in log space, so have to change expressions
accordingly. Multiplications get converted to additions in log space.  But to avoid numerical issues, *additions* of
probabilities are best converted by using the `log.sum.exp()` function (that should already be loaded), rather than
using `log( sum( exp()))` directly.

```R
# a generally useful function:
echo <- function( message, ... ) {
    cat( sprintf( message, ... ))
}


forward.backward <- function(
    emissions,              # matrix of emission log-probabilities (LxD)
    transitions,	   # (L-1)xDxD array of transition probabilities into each bin
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
```

## Computing transition probabilities

To define transition probabilities, we'll assume a model of exponential copy number tract length with expected length
controlled by a rate parameter `lambda`.  

We also pass in a list of between-site distances, but in our example all bins are the same size, so we will just pass
in a list of 1s and interpret parameters accordingly.  In other words, `lambda=1` would correspond to an expected tract length of 1600bp.

```R
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
```

## The Copy number HMM

The next function puts this all together to implement the CNV-finding HMM.
It takes:
* a matrix of per-site, per-sample copy number logliklihoods to be used as emission probabilities.
* prior probabilities for each copy number state
* and a lambda used to determine expected length of transition

It also takes a last parameter (`site.multipliers`) that can used to control for per-site variation, but we're not using it here.

```R
# Function cnv.hmm()
# Implement one iteration of the CNV-finding HMM for a matrix of N samples and L sites (bins)
#
# Return value is a list with several members
# - `prior` is the copy number prior, as passed in
# - `marginal.log.probabilities` are the the log-posterior probabilities of each copy number given the HMM (i.e. gamma in the HMM)
# - `expected.copy.numbers` are the expected copy numbers at each site, given the posteriors.
# - `emission.log.probabilities` the matrix of copy number lls, as passed in
# The total log-probability across all samples (this is the total log-likelihood of our model)
cnv.hmm <- function(
    copy.number.lls,
    prior,
    lambda,
    site.multipliers
) {

    # sanity checks
    stopifnot( dim( copy.number.lls )[3] == length( prior ))

    # Get a list of samples
    samples = dimnames(copy.number.lls)[[2]]
    N = length( samples )
    L = nrow( copy.number.lls) # number of sites

    echo( "Running HMM for %d samples at %d sites...\n", N, L )
    
    # prior is assumed to be for copy number states 0 ... K
    copy.numbers = 0:(length(prior)-1)
    log.prior = log( prior )

    result = list(
        prior = prior,
        # This array reports the posterior probability of each
        # copy number state at each site for each sample
        # under the HMM model
        marginal.log.probabilities = array(
            NA,
            dim = c( N, L, length( copy.numbers )),
            dimnames = list(
                samples,
                rownames( copy.number.lls ),
                sprintf( "cn=%d", copy.numbers )
            )
        ),
        expected.copy.numbers = array(
            NA,
            dim = c( L, N ),
            dimnames = list(
                rownames(copy.number.lls ),
                samples
            )
        ),
        emission.log.probabilities = copy.number.lls,
        total.log.probability = 0
    )
    transitions = compute.transitions( rep( 1, L-1 ), lambda, log.prior )
    for( i in 1:N ) {
        emissions = copy.number.lls[,i,]
        fb = forward.backward( emissions, transitions, log.prior )
        result$total.log.probability = result$total.log.probability + log.sum.exp( fb$alpha[1,] + fb$beta[1,] )
        result$marginal.log.probabilities[i,,] = fb$gamma
        result$expected.copy.numbers[,i] = exp(fb$gamma) %*% copy.numbers
        echo( "." )
        if( i %% 50 == 0 ) {
            echo( "\n" )
        }
    }
    echo( "\nok\n" )
    return( result ) ;
}
```

## Run the HMM

Ok we are ready to run!  Let's set it going:

```R
result = cnv.hmm(
    # Copy number log-likelihoods computed earlier.
    copy.number.lls,

    # Prior on each copy number state
    # States are assumed to start at zero (no copies) and go up to K
    # Here we put 90% prior on diploid state
    prior = c( 0.02, 0.02, 0.9, 0.02, 0.02, 0.02 ),

    # lambda = switch rate per bin
    # Each bin is 1600bp, so 1/10 expects 1 switch every 16kb and so on
    lambda = 1/20,

    # site multipliers affect values across samples for each site
    # these could be used to handle variation in coverage e.g. due to site-specific
    # mapping, sequence data rates, GC content etc.
    # Here we set these all to 1 for an initial run
    site.multipliers = rep( 1, nrow( X ))
)

plot.copy.numbers( result$expected.copy.number, title = "Expected copy number (HMM model)" )
```

Compare this with the previous plot based only on per-site coverage values for each sample - the new plot is much cleaner.  Also let's cluster samples again:

```R
o = hclust( # hierarchical clustering
	    dist( t(result$expected.copy.number) ) # of Euclidean distance matrix between samples
)$order
plot.copy.numbers( result$expected.copy.number[,o], title = "Expected copy number (HMM model, clustered)" )
```

## Further directions

Our model is getting better but there are lots of things we could do to improve this!

First, there's variation between sites e.g. sequence coverage variation due to genome GC content or mapping performance.
We could try to control for this by iteravely fitting per-site multipliers in the above (maybe using MCMC or another algorithm).

Second, having taken a first call of what the CNVs are, we now have more information to estimat ethe Gaussian parameters!  So we could go back and re-estimate, and iterate.

There are lots of other ways this model could be improved - challenge is to think of some.

