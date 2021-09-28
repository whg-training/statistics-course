library( numDeriv )

echo <- function( message, ... ) {
    cat( sprintf( message, ... ))
}

# Code to inspect a log-likelihood

# We've implemented functions with params as named lists
# To make the code below work I pass between this and a linearised
# representation of them in a parameter vector
destructure.params <- function( params ) {
    lengths = sapply( params, length )
    return(
        list(
            structure = list(
                names = names( params ),
                lengths = lengths
            ),
            parameter.vector = unlist( params )
        )
    )
}

restructure.params <- function( parameter.vector, structure ) {
    result = list()
    offset = 0
    for( i in 1:length( structure$names )) {
        name = structure$names[i]
        length = structure$lengths[[name]]
        result[[name]] = as.numeric( parameter.vector[offset+(1:structure$lengths[[name]]) ] )
        offset = offset + length
    }
    return( result ) 
}

find.mle <- function(
    ll,                     # loglikelihood function to optimise
    data,                   # data to pass as 1st argument
    params,                 # params to pass as 2nd argument (named list)
    bounds = NULL,          # any bounds
    ...                     # any extra data to be passed to ll()
) {
    destructured = destructure.params( params )
    echo( "Parameter structure is:\n" )
    print( destructured )
    
    dot.args = list(...)

    if( is.null( bounds )) {
        fn = function(p) {
            params = restructure.params( p, destructured$structure )
            if( "X" %in% names( dot.args )) {
                return( sum(ll( data, params, X = dot.args$X )))
            } else {
                return( sum(ll( data, params )))
            }
        }
    } else {
        fn = function(p) {
            params = restructure.params( p, destructured$structure )
            # To make optim() respect the bounds, we add a continuous penalty outside the
            # bounds, and clamp the params.
            penalty = 0
            for( name in names(params) ) {
                if( name %in% names( bounds )) {
                    bound = bounds[[name]]
                    if( !is.list( bound )) {
                        bound = list(bound)
                    }
                    for( i in 1:length( params[[name]] )) {
                        if( params[[name]][i] < bound[[i]][1] ) {
                            params[[name]][i] = bound[[i]][1] ;
                            penalty = penalty + (bound[[i]][1] - params[[name]][i])
                        } else if( params[[name]][i] > bound[[i]][2] ) {
                            params[[name]][i] = bound[[i]][2] ;
                            penalty = penalty + (params[[name]][i] - bound[[i]][2])
                        }
                    }
                }
            }
            if( "X" %in% names( dot.args )) {
                return( sum(ll( data, params, X = dot.args$X )) - penalty )
            } else {
                return( sum(ll( data, params )) - penalty )
            }
        }
    }
    D = sum( destructured$structure$lengths )
    if( D > 1 ) {
        method = "Nelder-Mead"
    } else {
        method = "Brent"
    }
    method = "Nelder-Mead"
    fit = optim(
        par = destructured$parameter.vector,
        fn = fn,
        control = list(
            trace = FALSE,
            fnscale = -1, # optim() maximises by default, so we invert here.
            warn.1d.NelderMead = FALSE
        ),
        hessian = FALSE,
#        method = "BFGS"
#        method = "L-BFGS-B"
        method = method
    )
    stopifnot( fit$converged == 0 )
    return(
        list(
            converged = (fit$converged == 0 ),
            mle = fit$par,
            hessian = hessian( fn, fit$par ),
            structure = destructured$structure
        )
    )
}

inspect.ll <- function(
    ll,
    data,
    params,
    selection = 1:length( params),
    bounds = NULL,
    asymptotic.approximation = TRUE,
    ... # Any extra data to be passed to ll() - currently must be called X
) {
    fit = find.mle( ll, data, params, bounds = bounds, ... )
    H = fit$hessian
    parameter.structure = fit$structure
    mle = fit$mle

    dot.args = list(...)

#    cat( "mle:\n" )
#    print(mle)
#    cat( "Hessian:\n" )
#    print(H)
#    cat( "Information:\n" )
    I = matrix( NA, nrow = length( mle ), ncol = length( mle ))
    tryCatch({
        I = solve(-H)
    }, error = function(e) {
       echo( "!! Error: hessian is singular.\n" )
       I = matrix( NA, length( mle )) 
    })
#    print(I)

    asymptotic.approximation = asymptotic.approximation && !is.na(sum(I))

    # insist on named params
    stopifnot( !is.null( names( params )))
    
    ll.fn = function(p) {
        params = restructure.params( p, parameter.structure )
        if( "X" %in% names( dot.args )) {
            return( sum(ll( data, params, X = dot.args$X )))
        } else {
            return( sum(ll( data, params )))
        }
    }
    
    D = length(mle)
    layout(
        matrix( 1:(2*D), byrow = T, ncol = 2 )
    )
    plot.data = data.frame()
    offset = 0
    for( name in parameter.structure$names ) {
        for( k in 1:parameter.structure$lengths[[name]] ) {
            echo( "Plotting %s[%d]...\n", name, k )
            i = offset + k
            se = sqrt(I[i,i])
            if( is.null( bounds ) || !name %in% names(bounds)) {
                if( !is.na( se )) {
                    size = min( 3*se, 100 )
                    parameter.bounds = c( mle[i] - size, mle[i] + size )
                } else {
                    parameter.bounds = c( mle[i] * 0.25, mle[i] * 4 )
                }
            } else {
                bound = bounds[[name]]
                if( is.list( bound )) {
                    parameter.bounds = bound[[k]]
                } else {
                    parameter.bounds = bound
                }
            }
        
            xvalues = seq( from = parameter.bounds[1], to = parameter.bounds[2], length = 100 )
            fn.values = sapply( xvalues, function(x) {
                p = mle;
                p[i] = x;
                return( ll.fn(p));
            })
            plot(
                xvalues, exp(fn.values),
                type = 'l',
                xlab = sprintf( "Parameter: %s [%d]", name, k ),
                ylab = "Likelihood (product over data)",
                bty = 'n',
                main = sprintf( "%s[%d]", name, k )
            )
            text(
                mle[i],
                exp(ll.fn( mle )),
                sprintf( "%s = %.3f, se = %.3f", name, mle[i], se ),
                pos = 3,
                cex = 0.75,
                font = 3,
                xpd = NA
            )
            if( asymptotic.approximation ) {
                # We want the approximation of parameter i by the conditional distribution
                # holding all other params at their mle values.
                # It is given by the formulae at https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
                conditional.mean = mle[i] # conditioning on mle does not affect maximum
                if( D > 1 ) {
                    Sigma11 = I[i,i,drop=F]
                    Sigma12 = I[i,-i,drop=F]
                    Sigma22 = I[-i,-i,drop=F]
                    Sigma22.inverse = solve(Sigma22)
                    conditional.variance = Sigma11 - Sigma12 %*% Sigma22.inverse %*% t(Sigma12)
                    conditional.se = sqrt( conditional.variance )
                } else {
                    conditional.mean = mle[i]
                    conditional.se = sqrt( I[i,i] )
                }
                gaussian.approximation = (
                    dnorm( xvalues, mean = conditional.mean, sd = conditional.se, log = T )
                    - dnorm( conditional.mean, mean = conditional.mean, sd = conditional.se, log = T )
                    + ll.fn( mle )
                )
                points( xvalues, exp(gaussian.approximation), type = 'l', lty = 2, col = 'gold3' )
            }
            grid()
            abline( v = mle[i], col = 'red', lty = 2 )

            if( asymptotic.approximation ) {
                legend(
                    "topleft",
                    legend = c( "ll", "approx." ),
                    lty = c( 1, 2 ),
                    col = c( "black", "gold3" ),
                    bty = 'n'
                )
            } else {
                legend(
                    "topleft",
                    legend = c( "ll" ),
                    lty = c( 1 ),
                    col = c( "black" ),
                    bty = 'n'
                )
            }
                    
            plot(
                xvalues, fn.values,
                type = 'l',
                xlab = sprintf( "Parameter: %s [%d]", name, k ),
                ylab = "Log-likelihood (summed over data)",
                bty = 'n',
                main = sprintf( "%s[%d] (log space)", name, k )
            )
            text(
                mle[i],
                ll.fn( mle ),
                sprintf( "%s = %.3f, se = %.3f", name, mle[i], se ),
                pos = 3,
                cex = 0.75,
                font = 3,
                xpd = NA
            )
            if( asymptotic.approximation ) {
                points( xvalues, gaussian.approximation, type = 'l', lty = 2, col = 'gold3' )
            }

            grid()
            abline( v = mle[i], col = 'red', lty = 2 )
        }
        offset = offset + parameter.structure$lengths[[name]]
    }
    return( list(
        H = H,
        I = I,
        mle = restructure.params( mle, parameter.structure ),
        standard.errors = sqrt( diag( I ))
    ))
}

restrict <- function( ll, to ) {
    result = function( data, params, ... ) {
        params[names(to)] = to
        ll( data, params, ... )
    }
    return( result )
}

# EXAMPLES

# Gaussian example, 2 observations
inspect.ll(
    gaussian.ll,
    data = c( 1, 2 ),                       # 2 observations
    params = list( mean = 1, sigma2 = 0.2 ) # starting values
)

# Gaussian with fixed variance
inspect.ll(
    restrict( gaussian.ll, to = list( sigma2 = 1 )),
    data = c( 1, 2 ),
    params = list( mean = 1 )
)

# 2x2 table example
source("binomial.ll.R" )
source("table.ll.R" )
table.ll( TBL, params = list(
    theta = c( 0.5, 0.5 )
))
A = matrix( c( 2, 5, 7, 8 ), byrow = T, ncol = 2 )
inspect.ll(
    table.ll,
    A,
    params = list(
        theta = c( 0.5, 0.5 )
    )
)

A = matrix( c( 1, 100, 0, 10 ), byrow = T, nrow = 2 )
inspect.ll(
    reparameterised.table.ll,
    A,
    params = list(
        mu = 0.5, log.or = 0.5
    )
)


inspect.ll(
    reparameterised.table.ll,
    TBL,
    params = list(
        theta = 0.5,
        log.or = 0
    )
)

inspect.ll(
    restrict( binomial.ll, to = list( n = 1 )),
    c( 1, 1, 0 ),
    params = list( p = 0.2 ),
    bounds = list( p = c( 0, 1 ))
)

inspect.ll(
    restrict( binomial.ll, to = list( n = 10 )),
    rep(1,100),
    params = list( p = 0.2 )
)

inspect.ll(
    restrict( binomial.ll, to = list( n = 10 )),
    rep(1,100),
    params = list( p = 0.2 ),
    bounds = list( p = c( 0, 0.4 ) )
)
