# Requires binomial.ll first!
table.ll <- function( data, params = list( theta = c( 0.5, 0.5 )) ) {
    # For a 2x2 table we assume binomial sampling in rows
    # So sum the lls over the rows
    return (
        binomial.ll(
            x = data[1,2],
            params = list(
                n = rowSums( data )[1],
                p = params$theta[1]
            )
        )
        + binomial.ll(
            data[2,2],
            params = list(
                n = sum( data[2,] ),
                p = params$theta[2]
            )
        )
    )
}

logistic <- function( x ) {
	exp(x) / ( 1 + exp(x) )
}

reparameterised.table.ll <- function( data, params = list( theta = 0.5, log.or = 0 )) {
	# params[1] is log-odds of baseline frequency
	# params[2] is log odds ratio
	theta = c(
		logistic( params$theta ),
		logistic( params$theta + params$log.or )
	)
	return( table.ll( data, params = list( theta = theta )))
}

