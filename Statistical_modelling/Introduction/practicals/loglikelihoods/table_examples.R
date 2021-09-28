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
