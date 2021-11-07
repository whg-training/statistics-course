library( tidyverse )
data = read_csv( "../o_bld_grp/o_bld_grp.csv" )
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

svg( filename = "solutions/Tanzania_o_blood_group_inference.svg", width = 4, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_unnormalised_posterior( 41, 102 )
dev.off()

layout( matrix( 1:8, ncol = 2 ))
countries = unique( data$country )
for( country in countries ) {
	par( mar = c(4,2,1,0.1) )
	w = which( data$country == country )
	counts = table( data$O.bld.grp[w] )
	plot_unnormalised_posterior( counts[1], counts[1] + counts[2] )
	legend( "topleft", legend = c( country, sprintf( "(%d, %d)", counts[1], counts[2] )), bty = 'n' )
}

layout( matrix( 1:24, ncol = 3 ))
ethnicities = unique( data$ethnicity )
for( ethnicity in ethnicities ) {
	par( mar = c(2,2,1,0.1) )
	w = which( data$ethnicity == ethnicity )
	counts = table( data$O.bld.grp[w] )
	plot_unnormalised_posterior( counts[1], counts[1] + counts[2] )
	legend( "topleft", legend = c( ethnicity, sprintf( "(%d, %d)", counts[1], counts[2] )), bty = 'n' )
}

plot_normalised_posterior <- function( k, n ) {
	f <- function( y ) { return( binomial.likelihood( k, n, y ) ) ; }
	denominator = integrate( f, 0, 1 )$value

	x = seq( from = 0, to = 1, by = 0.01 )
	plot(
		x, binomial.likelihood( k = k, n = n, theta = x ) / denominator,
		type = 'l',
		xlab = "O blood group frequency (θ)",
		ylab = "Normalised posterior",
		yaxt = 'n',
		bty = 'n'
	)
	grid()
	abline( v = k/n, col = 'red' )
	axis( 2, las = 1 )
}

svg( filename = "solutions/Tanzania_o_blood_group_normalised.svg", width = 4, height = 2 )
par( mar = c(4,5,1,0.1) )
plot_normalised_posterior( 41, 102 )
dev.off()

