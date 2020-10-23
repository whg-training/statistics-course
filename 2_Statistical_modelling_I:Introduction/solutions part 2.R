################################
# Functions from last week
logistic <- function( x ) {
	exp(x) / ( 1 + exp(x) )
}

binomial.ll <- function( y, n, p ) {
  result = lchoose( n, y ) + y*log(p) + (n-y) * log(1-p) ;
  return( result )
}

find.mle <- function( outcome, design.matrix, loglikelihood ) {
	fit = optim(
		par = c( 0, 0 ),
		fn = function( params ) {
			loglikelihood( outcome, design.matrix, params )  
		},
		# tell optim() to maximise, not minimise.
		control = list(
			fnscale = -1,
			trace = TRUE
		)
	)
	if( fit$convergence != 0 ) {
		return( NA )
	} else {
		return( fit$par ) # estimate
	}
}
""

# Challenge #1
logistic.regression.ll <- function( Y, design.matrix, params ) {
	# outcome is an Nx1 vector of 1's (e.g. cases) and 0's (e.g. controls)
	# design.matrix is an Nxd matrix of predictors and covariates
	# params is a dx1 column vector of parameters
	predictors = design.matrix %*% params # Nb. %*% means matrix multiplication
	lls = binomial.ll( Y, 1, logistic( predictors ))
	return( sum( lls, na.rm = T ))
}

# Compute the 2nd derivative of logistic regression ll
# Implemented by GB - not part of the practical.
logistic.regression.ddll <- function( Y, design.matrix, params ) {
	# should be one 
	predictors = design.matrix %*% params # matrix multiplication
	probs = logistic( predictors )
	f = probs * ( 1 - probs )
	result = 0.0
	for( i in 1:nrow(design.matrix)) {
		if( !is.na( sum( design.matrix[i,]) )) {
			result = result - kronecker( t( design.matrix[i,,drop = F] ), design.matrix[i,,drop=F] ) * f[i]
		}
	}
	return( result )
}

plot.loglikelihood <- function( Y, design.matrix, loglikelihood ) {
	# loglikelihood is a log-likelihood function
	# coefficients is the result of summary(fit)$coeffs where fit = glm( ... )
	# Get the number of parameters

	# mean-centre columns.  This is to allow the normal approximation
	# used below (based on marginal standard errors) to look right
	for( i in 2:ncol( design.matrix )) {
		design.matrix[,i] = design.matrix[,i] - mean( design.matrix[,i], na.rm = T )
	}
	
	fit = glm( Y ~ design.matrix - 1, family = "binomial" )
	coefficients = summary(fit)$coeff
	d = nrow( coefficients )
	betas = coefficients[,'Estimate']
	max.ll = loglikelihood( Y, design.matrix, betas )
	par( mfrow = c( (d+1)/2, 2 ), mar = c( 2, 3, 2, 1 ))
	for( i in 1:d ) {
		beta_hat = betas[i]
		se = coefficients[i,'Std. Error']
		at = seq( from = beta_hat - pmin( se, 2 ) * 5, to = beta_hat + pmin( se, 2 ) * 5, by = 0.01 )
		evaluations = sapply( at, function(x) {
			params = betas
			params[i] = x
			return( loglikelihood( Y, design.matrix, params ) )
		})
		plot(
			at,
			evaluations,
			type = 'l',
			main = gsub( "design.matrix", "", rownames( coefficients )[i], fixed = T ),
			xlab = "",
			ylab = "ll",
			lwd = 2
		)
		grid()
		
		# Overlay a gaussian to see what it looks like
		gaussian.evaluations = -0.5 * ( at - beta_hat )^2 / se^2
		max.gaussian = 0
		points(
			at,
			gaussian.evaluations - max.gaussian + max.ll,
			type = 'l',
			lty = 3,
			col = "red"
		)
		abline( v = beta_hat, col = 'red' )
		text(
			beta_hat, min( evaluations ), pos = 4,
			sprintf(
				"%.2f (%.2f - %.2f)",
				beta_hat,
				beta_hat - 1.96 * se,
				beta_hat + 1.96 * se
			),
			col = "red"
		)
	}
}

fixed.effect.meta <- function(
	betas,
	ses
) {
	W = 1 / sum( 1/ses^2 )
	scaled_betas = betas / ses^2
	return( list(
		meta.beta = W * sum( scaled_betas ),
		meta.se = sqrt(W)
	))
}

