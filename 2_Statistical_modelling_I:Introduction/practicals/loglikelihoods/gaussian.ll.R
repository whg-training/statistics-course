gaussian.ll <- function( x, params = list( mean = 1, sigma2 = 1 )) {
    mean = params$mean
    sigma2 = params$sigma2

    if( min( params$sigma2 ) < 0 ) {
        return( -Inf )
    }

    result = -0.5 * log( 2*pi*sigma2 ) - (x-mean)^2/(2*sigma2)
    return( result )
}
