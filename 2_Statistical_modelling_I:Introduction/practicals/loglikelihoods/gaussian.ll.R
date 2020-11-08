gaussian.ll <- function( x, params = list( mean = 1, sigma2 = 1 )) {
    # Handle negative sigmas
    if( min( params$sigma2 ) < 0 ) {
        return( -Inf )
    }
    with( params, -0.5 * ( log(sigma2) + log( 2*pi ) + (x-mean)^2/sigma2 ) )
    # or use built-in function:
    #dnorm( x, mean = 1, sd = sqrt( params$variance ), log = T )
}
