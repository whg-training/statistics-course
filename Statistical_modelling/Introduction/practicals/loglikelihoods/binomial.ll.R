# My implementation of binomial log-likelihood
# Note this allows a continuous value of n!
# (The maths does not disallow this, and it makes optimisation easier, so why not?)
#
# NoteL: like dbinom this allows multiple x values, and also multiple n and p values.
# n and p parameters must have the same length, which must either be length 1
# or the same length as x
# 
binomial.ll <- function(
    x,                                      # number of successes
    params = list(
        n = 1,                              # number of trials
        p = 0.5                             # success probability
    )
) {
    # Either:
    # n and p have length 1
    # OR:
    # n and p have same length as x.
    # otherwise we quit.
    n = params$n
    p = params$p
    if( length(n) == 1 ) {
        n = rep( n, length(x) )
    }
    if( length(p) == 1 ) {
        p = rep( p, length(x) )
    }
    stopifnot( length(n) == length(p) )
    stopifnot( length(x) == length(n) )

    # Binomial distribution is:
    # (n choose x ) * p^x * (1-p)^(n-x)
    # Log-likelihood is:
    # log(n choose x) + x * log(p) + (n-x) * log(1-p)
    # Any 0^0 s in the above should be interpreted as 1 i.e. log = 0
    # Therefore we avoid adding the terms when x=0 and x=n
    result = lchoose( n, x ) ;
    result[ x>0 ] = result[x>0] + x[x>0] * log(p[x>0])
    result[ x<n ] = result[x<n] + (n[x<n] - x[x<n]) * log(1-p[x<n])

    # Handle degenerate cases
    result[ p < 0 | p > 1 ] = -Inf
    result[ n < 0 ] = -Inf
    result[ x<0 ] = -Inf
    result[ x>n ] = -Inf
    return(result)
}

# Here's a version of binomial.ll that uses dbinom to do the computation
# It is simpler, but rounds n and this makes optimisation over n harder
binomial.ll.dbinom_version <- function(
    x,
    params = list( n = 1, p = 0.5 )
) {
    dbinom( x, size = round(params$n), prob = params$p, log = T )
}

