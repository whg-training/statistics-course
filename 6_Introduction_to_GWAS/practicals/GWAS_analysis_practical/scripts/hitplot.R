# Get a region

hitplot <- function( focus.snp, data, genes, ld, margin = 200000 ) {
# Get a region around the focus SNP
  focus.pos = data$BP[ data$SNP == focus.snp ]
  region = c( focus.pos - margin, focus.pos + margin )

  data = data[ data$BP >= region[1] & data$BP <= region[2], ]

  # Use heat.colors() to colour points from red to white
  # according to LD
  region.ld = ld[ ld$SNP_A == focus.snp & ld$BP_B >= region[1] & ld$BP_B <= region[2] & !is.na( ld$R2_bin ), ]
  data$colour = "black"
  data$colour[ match( region.ld$SNP_B, data$SNP ) ] = heat.colors(10)[ as.integer( region.ld$R2_bin ) ]

  # Make all points round
  if( !'shape' %in% names( data ) ) {
	  data$shape = 20
  }

  # Get a region around the focus SNP
  layout.matrix = matrix(c(1,2,3), ncol=1 )
  layout( layout.matrix, heights = c(1,0.4,0.4) )

  # Plot the association signal
  par( mar = c( 0.1,4,1,1 ))
  plot(
	  data$BP, -log10( data$P ),
	  col = data$colour, pch = data$shape,
	  xlim = region,
	  xlab = '', ylab = "-log10( P-value )",
	  xaxt = 'n',
	  mar=c(0.1,4,1,1)
  )

  # Since this is a stacked plot, it's nice to plot the x axis at
  # the top and bottom of the whole plot.
  # To do that we get the x axis tick points and manually place
  # the axis using the axis() command.  (The tcl and mgp arguments
  # are there to make subtle adjustments to the positioning of the ticks and labels).
  ticks = axTicks(3)
  axis( side = 3, at = ticks, labels = sprintf( "%.1fMb", as.numeric( ticks ) / 1000000 ), tcl = 0.5, mgp = c( 3,0.1,0 ) )

  # Add a legend
  legend(
	  "topright",
	  col = c( "black", heat.colors(10) ),
	  pch = 20,
	  legend = c( "r²<0.1", sprintf( "r²≥%.1f", seq( from = 0.1, 0.9, by = 0.1 ) ) )
  )
  grid()

  # Plot genes.  We remove the border with bty='n'
  par( mar = c( 0.1,4,0.1,1 ))
  plot.genes( genes, region, bty = 'n', yaxt = 'n', xaxt = 'n' )
  grid( ny = NA )

  # Plot the recombination map
  par( mar = c( 1.5,4,0.1,1 ))
  plot( genetic_map$position, genetic_map$ COMBINED_rate.cM.Mb., type = 'l', xlim = region, xaxt = 'n', ylab = "cM/Mb", mar = c( 1,1,0.1,2 ) )
  # Plot the x axis again
  axis( side = 1, at = ticks, labels = sprintf( "%.1fMb", as.numeric( ticks ) / 1000000 ), tcl = 0.5, mgp = c( 3,0.5,0 ) )
  grid()
}
