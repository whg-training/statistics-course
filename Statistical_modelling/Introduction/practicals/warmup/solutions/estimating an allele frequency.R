library( tidyverse )
data = read_csv( "o_bld_grp.csv" )
table( data$country, data$O.bld.grp )



binomial.likelihood <- function( k, n, theta ) {
    return( dbinom( x = k, size = n, prob = theta ) )
}

plot_unnormalised_posterior <- function( k, n ) {
	x = seq( from = 0, to = 1, by = 0.01 )
	plot(
		x, binomial.likelihood( k = k, n = n, theta = x ),
		type = 'l',
		xlab = "O blood group frequency (θ)",
		ylab = "Unnormalised\nposterior",
		yaxt = 'n',
		bty = 'n'
	)
	grid()
	abline( v = k/n, col = 'red' )
	axis( 2, las = 1 )
}

counts = table( data$country, data$O.bld.grp )['Tanzania',]
k = counts[2]
n = sum(counts)

svg( filename = "solutions/Tanzania_o_blood_group_inference.svg", width = 4, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_unnormalised_posterior( k, n )
dev.off()

plot.countries <- function( data ) {
	countries = unique( data$country )
	layout( matrix( 1:(ceiling(length(countries)/3)*3), ncol = 3 ))
	for( country in countries ) {
		par( mar = c(4,2,1,0.1) )
		w = which( data$country == country )
		counts = table( data$O.bld.grp[w] )
		k = counts[2]
		n = sum(counts)
		plot_unnormalised_posterior( k, n )
		legend( "topleft", legend = c( country, sprintf( "(%d, %d)", counts[1], counts[2] )), bty = 'n' )
	}
}

plot.ethnicities <- function( data ) {
	ethnicities = unique( data$ethnicity )
	layout( matrix( 1:(ceiling(length(ethnicities)/3)*3), ncol = 3 ))
	for( ethnicity in ethnicities ) {
		par( mar = c(2,2,1,0.1) )
		w = which( data$ethnicity == ethnicity )
		counts = table( data$O.bld.grp[w] )
		plot_unnormalised_posterior( counts[2], counts[1] + counts[2] )
		legend( "topleft", legend = c( ethnicity, sprintf( "(%d, %d)", counts[1], counts[2] )), bty = 'n' )
	}
}

plot_normalised_posterior <- function( k, n ) {
	f <- function( y ) { return( binomial.likelihood( k, n, y ) ) ; }
	denominator = integrate( f, 0, 1 )$value

	x = seq( from = 0, to = 1, by = 0.01 )
	plot(
		x, binomial.likelihood( k = k, n = n, theta = x ) / denominator,
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

svg( filename = "solutions/Tanzania_o_blood_group_posterior.svg", width = 4, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_normalised_posterior( k, sum(counts) )
dev.off()

svg( filename = "solutions/Tanzania_o_blood_group_posterior+beta.svg", width = 4.5, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_normalised_posterior( k, sum(counts) )
x = seq( from = 0, to = 1, by = 0.01 )
points( x, dbeta( x, shape1 = k+1, shape2 = n-k+1 ), type = 'l', lty = 2, col = 'blue', lwd = 2 )
legend( "topleft", col = c( "black", "blue" ), lty = c( 1, 2 ), , lwd = c( 1, 2 ), legend = c( "Numerical", sprintf( "Beta(%d,%d)", counts[1]+1, counts[2]+1 )), bty = 'n' )
dev.off()


plot_beta_density <- function( k, n ) {
	x = seq( from = 0, to = 1, by = 0.01 )
	plot( x, dbeta( x, shape1 = k+1, shape2 = (n-k)+1 ), type = 'l', ylab = "Beta density", xlab = "O blood group frequency (θ)" )
	grid()
	abline( v = k/n, col = 'red' )
}

svg( filename = "solutions/Tanzania_o_blood_group_beta_posterior.svg", width = 4, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_beta_density( k, sum(counts) )
dev.off()

compute.credible.interval <- function( k, n, mass = 0.95 ) {
	tail = 1 - mass
	return( qbeta( c( lower = tail/2, median = 0.5, upper = 1.0-(tail/2) ), shape1 = k+1, shape2 = n-k+1 ))
}

# table() doesn't work well if one of the counts is zero
compute.counts <- function( data ) {
	return( c(
		nonO = length( which( data$O.bld.grp == 0 )),
		O = length( which( data$O.bld.grp == 1 ))
	))
}
summarise <- function( data.subset, name ) {
	counts = compute.counts( data.subset )
	credible = compute.credible.interval( counts[2], sum(counts ))
	return( c(
		list( name = name, nonO = counts[1], O = counts[2], estimate = counts[2]/sum(counts) ),
		credible
	))
}

# Using purr's map functions.
map_dfr( countries, function( country ) { summarise( data[ data$country == country, ], country ) } )
map_dfr( ethnicities, function( ethnicity ) { summarise( data[ data$ethnicity == ethnicity, ], ethnicity ) } )

plot.ethnicities <- function( data, prior.counts = c( 0, 0 ) ) {
	ethnicities = unique( data$country_ethnicity )
	layout( matrix( 1:(ceiling(length(ethnicities)/3)*3), ncol = 3 ))
	for( ethnicity in ethnicities ) {
		par( mar = c(2,2,1,0.1) )
		w = which( data$country_ethnicity == ethnicity )
		summary = summarise( data[w,], ethnicity )
		counts = compute.counts( data[w,] )
		plot_unnormalised_posterior( counts[2], counts[1] + counts[2] + sum(prior.counts) )
		rect( xleft = summary$lower, xright = summary$upper, ybottom = 0, ytop = 1, col = rgb(0,0,0,0.1), border = NA )
		legend( "topleft", legend = c( ethnicity, sprintf( "(%d, %d)", counts[1], counts[2] )), bty = 'n' )
	}
}

svg( filename = "solutions/all_ethnicities_o_blood_group_beta_posterior.svg", width = 8, height = 8 )
par( mar = c(4,5,1,0.1) )
plot.ethnicities( data )
dev.off()
