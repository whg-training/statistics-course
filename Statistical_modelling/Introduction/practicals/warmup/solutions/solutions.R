# table() doesn't work well if one of the counts is zero
# so write compute.counts instead
compute.counts <- function( data ) {
	return( c(
		nonO = length( which( data$O.bld.grp == 0 )),
		O = length( which( data$O.bld.grp == 1 ))
	))
}
compute.credible.interval <- function( k, n, mass = 0.95 ) {
	tail = 1 - mass
	return( qbeta( c( lower = tail/2, median = 0.5, upper = 1.0-(tail/2) ), shape1 = k+1, shape2 = n-k+1 ))
}

summarise <- function( data.subset, name, prior.counts = c( 0, 0 )) {
	counts = compute.counts( data.subset )
	credible = compute.credible.interval( counts[2] + prior.counts[1], sum(counts) + sum(prior.counts))
	return( c(
		list( name = name, nonO = counts[1], O = counts[2], estimate = counts[2]/sum(counts) ),
		credible
	))
}
plot_normalised_posterior <- function( k, n ) {
	x = seq( from = 0, to = 1, by = 0.01 )
	plot(
		x, dbeta( x, shape1 = k+1, shape2 = n-k+1),
		type = 'l',
		xlab = "O blood group frequency (θ)",
		ylab = "posterior",
		yaxt = 'n',
		bty = 'n'
	)
	grid()
	abline( v = k/n, col = 'red' )
	axis( 2, las = 1 )
}
plot.bygroup <- function( data, column, prior.counts = c( 0, 0 )) {
	levels = unique( data[[column]] )
	print(levels)
	layout(
		matrix(
			1:(ceiling(length(levels)/3)*3),
			ncol = 3
		)
	)
	for( level in levels ) {
		par( mar = c(2,2,1,0.1) )
		w = which( data[,column] == level )
		summary = summarise( data[w,], level, prior.counts )
		plot_normalised_posterior( summary$O + prior.counts[2], summary$O + summary$nonO + sum(prior.counts))
		rect( xleft = summary$lower, xright = summary$upper, ybottom = 0, ytop = 1, col = rgb(0,0,0,0.1), border = NA )
		legend( "topleft", legend = c( level, sprintf( "(%d, %d)", summary$nonO, summary$O )), bty = 'n' )
	}
}

