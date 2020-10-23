tables = list(
	band_et_al = matrix(
		c( 3420, 3233, 3925, 2738 ),
		byrow = T, ncol = 2,
		dimnames = list(
			c( "non-O", "O" ),
			c( "controls", "cases" )
		)
	),
	egan_et_al = matrix(
		c( 1965, 1, 707, 17 ),
		byrow = T, ncol = 2,
		dimnames = list(
			c( "G", "C" ),
			c( "unexposed", "exposed" )
		)
	)
)

binomial.likelihood <- function( y, n, p ) {
  result = choose( n, y ) * p^y * (1-p)^(n-y) ;
  return( result )
}

binomial.ll <- function( y, n, p ) {
  result = lchoose( n, y ) + y*log(p) + (n-y) * log(1-p) ;
  return( result )
}

table.ll <- function( data, params ) {
	return (
		binomial.ll( data[1,2], sum( data[1,] ), params[1] )
		+ binomial.ll( data[2,2], sum( data[2,] ), params[2] )
	)
}

logistic <- function( x ) {
	exp(x) / ( 1 + exp(x) )
}

reparameterised.table.ll <- function( data, params ) {
	# params[1] is frequency
	# params[2] is log Odds ratio
	thetas = c(
		logistic( params[1] ),
		logistic( params[1] + params[2] )
	)
	return( table.ll( data, thetas ))
}

find.mle <- function( data ) {
	fit = optim(
		par = c( 0, 0 ),
		fn = function( params ) {
			reparameterised.table.ll( data, params )  
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

plot.ll <- function( data, maximum.likelihood.estimate ) {
	# This solves Challenge #6
	mu_hat = maximum.likelihood.estimate[1]
	beta_hat = maximum.likelihood.estimate[2]
	se = sqrt( sum(1/(data+0.5) ))
	range = beta_hat + c( -50, 50 ) # may need to adjust
	at = seq( from = range[1], to = range[2], by = 0.01 )
	ll.evaluations = sapply(
		at,
		function(x) {
			reparameterised.table.ll(
				data,
				c( mu_hat, x )
			) }
	)

	# Plot likelihood on log scale
	plot(
		at,
		ll.evaluations,
		type = "l",
		xlab = "log odds Ratio", ylab = "log-likelihood",
		ylim = c( -10000, 0 )
	)
	max.ll = reparameterised.table.ll(
		data,
		maximum.likelihood.estimate
	) 
	max.normal.ll = dnorm( beta_hat, mean = beta_hat, sd = se, log = T )
	points(
		at,
		dnorm( at, mean = beta_hat, sd = se, log = TRUE ) - max.normal.ll + max.ll,
		type = 'l', lty = 3
	)
	
	grid()

	# annotate with the maximum likelihood estimate
	abline( v = beta_hat, col = "red" )
	text( beta_hat, -1000, "mle", col = "red", pos = 2 )

	# Plot likelihood on non-log scale
	plot(
		at,
		exp(ll.evaluations),
		xlim = beta_hat + c( -1, 1 ),  # zoom in a bit
		type = "l",
		xlab = "log odds Ratio", ylab = "likelihood"
	)

	points(
		at,
		exp( dnorm( at, mean = beta_hat, sd = se, log = TRUE ) - max.normal.ll + max.ll ),
		type = 'l',
		lty = 3
	)
	grid()

	# annotate with the maximum likelihood estimate
	abline( v = beta_hat, col = "red" )
}

simulate.tables <- function(
	log.or = -0.3,
	control.frequency = c( 0.5 ),
	sample.sizes = c( 10, 100, 1000, 1000, 5000 )
) {
	result = list()
	for( theta1 in control.frequency ) {
		for( N in sample.sizes ) {
			mu = log( theta1 / (1 - theta1 ))
			beta = log.or
			M = matrix( NA, nrow = 2, ncol = 2 )
			M[1,2] = rbinom( 1, N, prob = logistic( mu ) )
			M[2,2] = rbinom( 1, N, prob = logistic( mu + beta ) )
			M[1,1] = N - M[1,2]
			M[2,1] = N - M[1,2]
			result[[ sprintf( "theta1=%.2f/N=%d", theta1, N ) ]] = M
		}
	}
	return( result )
}

plot.simulated.tables <- function( simulated.data ) {
	par(
		mfrow = c( length(simulated.data), 2 ),
		mar = c( 2, 2, 1, 1 )
	)
	for( i in 1:length( simulated.data )) {
		plot.ll(
			simulated.data[[i]],
			find.mle( simulated.data[[i]] )
		)
	}
}

plot.simulated.tables( simulate.tables() )
