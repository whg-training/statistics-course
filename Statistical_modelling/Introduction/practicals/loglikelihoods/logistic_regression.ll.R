logistic <- function( x ) {
	exp(x) / ( 1 + exp(x) )
}

logistic_regression.ll <- function( Y, params = list( beta = c() ), X ) {
	# outcome is an Nx1 vector of 1's (e.g. cases) and 0's (e.g. controls)
	# design.matrix is an Nxd matrix of predictors and covariates
	# params is a dx1 column vector of parameters
	predictors = X %*% params$beta # Nb. %*% means matrix multiplication
	binomial.ll( Y, params = list( n = 1, p = logistic( predictors )))
}

# Compute the 2nd derivative of logistic regression ll
# Implemented by GB - not part of the practical.
logistic_regression.ddll <- function( Y, design.matrix, params ) {
	# should be one 
	predictors = design.matrix %*% params # matrix multiplication
	probs = logistic( predictors )
	f = probs * ( 1 - probs )
	result = 0.0
	for( i in 1:nrow(design.matrix)) {
		if( !is.na( sum( design.matrix[i,]) )) {
			result = result - kronecker( t( design.matrix[i,,drop = F] ), design.matrix[i,,drop=F]) * f[i]
		}
	}
	return( result )
}

logistic_regression.snptest.ll <- function( Y, params = list( beta = c() ), X ) {
	# Last three columns of X assumed to be 
	Z = X[, 1:(ncol(X)-2)]
	terms = matrix( NA, nrow = nrow(X), ncol = 3 )
	for( i in 1:3 ) {
		Z[,ncol(Z)] = X[,ncol(X)-3+i]
		predictors = Z %*% params$beta # Nb. %*% means matrix multiplication
		terms[,i] = binomial.ll( Y, params = list( n = 1, p = logistic( predictors )))
	}
	print( head( terms ))
	terms = terms + log( X[,(ncol(X)-2):ncol(X)])
	print( head( terms ))
	result = log( rowSums( exp(terms) ) )
	return( result )
}

