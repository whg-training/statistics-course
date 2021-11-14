simulate.population <- function(
  L,       # number of bi-allelic variants
  N,       # population size
  G,       # number of generations
  H = 10,  # starting number of haplotypes
  initial.frequency = 0.2 # initial variant frequency
) {
  result = list(
    L = L, H = H, N = N, G = G,
    haplotypes = NA,
    population = NA
  )
  haplotypes = matrix( NA, nrow = H, ncol = L )
  haplotypes[,] = rbinom( n = H*L, size = 1, prob = initial.frequency )
  population = matrix(
    NA, nrow = G, ncol = N,
    dimnames = list(
      sprintf( "generation=%d", 1:G ),
      sprintf( "individual=%d", 1:N )
    )
  )
  population[1,] = sample( 1:H, size = N, replace = TRUE )
  
  evolve <- function( currentGeneration ) {
    parents = sample( 1:N, size = N, replace = TRUE )
    return( currentGeneration[parents] )
  }
  
  for( generation in 2:G ) {
    population[generation,] = evolve( population[generation-1,] )
  }
  result$haplotypes = haplotypes
  result$population = population
  
  return(result)
}

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = "", ylab = "" ) {
  # this function plots a blank canvas
  plot(
    0, 0,
    col = 'white', # draw points white
    bty = 'n',     # no border
    xaxt = 'n',    # no x axis
    yaxt = 'n',    # no y axis
    xlab = xlab,   # no x axis label
    ylab = ylab,   # no x axis label
    xlim = xlim,
    ylim = ylim
  )
}

compute.ld <- function( haplotypes ) {
  cor( haplotypes )
}

compute.haplotype.frequencies <- function( sim ) {
  # We compute the frequency of each variant at each generation
  result = matrix(
    nrow = sim$G,
    ncol = sim$H,
    dimnames = list(
      sprintf( "generation=%d", 1:sim$G ),
      sprintf( "haplotype=%d", 1:sim$H )
    )
  )
  for( i in 1:sim$G ) {
    currentGeneration = sim$population[i,]
    result[i,] = sapply( 1:sim$H, function(h) { length( which( currentGeneration == h )) } ) / sim$N
  }
  return( result )
}

compute.variant.frequencies <- function( sim ) {
  # We compute the frequency of each variant at each generation
  result = matrix(
    nrow = sim$G,
    ncol = sim$L,
    dimnames = list(
      sprintf( "generation=%d", 1:sim$G ),
      sprintf( "variant=%d", 1:sim$L )
    )
  )
  for( i in 1:sim$G ) {
    currentGeneration = sim$haplotypes[sim$population[i,],]
    result[i,] = colSums( currentGeneration ) / sim$N
  }
  return( result )
}

plot.haplotypes <- function( haplotypes ) {
  image(
    t(haplotypes),
    x = 1:ncol(haplotypes),
    y = 1:nrow(haplotypes),
    xlab = "variant",
    ylab = "haplotype",
	breaks = c( -0.01, 0.5, 1 ),
	col = c( "#FFFFC8", "#7D0025" )
#	col = c( "#4B0055", "#FDE333" )
  )
}

plot.simulation <- function( sim ) {
  layout(
    matrix(
      c( 1, 1,
         2, 4,
         3, 5,
         6, 7
      ),
      byrow = T,
      nrow = 4
    ),
    heights = c( 1/10, 2/5, 2/5, 2/5 )
  )
  par( mar = c( 1, 3, 1, 1 ))

  # First let's write a legend
  blank.plot()
  legend(
    "center",
    legend = c(
      sprintf( "number of variants (L) = %d", sim$L ),
      sprintf( "number of haplotypes (H) = %d", sim$H ),
      sprintf( "number of individuals (N) = %d", sim$N ),
      sprintf( "number of generations (G) = %d", sim$G )
    ),
    bty = 'n',
    xpd = NA,
    ncol = 2
  )

  par( mar = c( 3, 3, 3, 1 ))
  # Plot starting and ending haplotypes
  plot.haplotypes( sim$haplotypes[sim$population[1,],] )
  mtext( "Variants", 1, line = 2, cex = 0.5 )
  mtext( "Haplotypes", 2, line = 2, cex = 0.5 )
  mtext(
    sprintf(
      "1st generation (%d segregating haplotypes)",
      length( unique( sim$population[1,] ))
    ),
    3,
    line = 1
  )
  plot.haplotypes( sim$haplotypes[sim$population[sim$G,],] )
  mtext( "Variants", 1, line = 2, cex = 0.5 )
  mtext( "Haplotypes", 2, line = 2, cex = 0.5 )
  mtext(
    sprintf(
      "Generation %d (%d segregating haplotypes)",
      sim$G,
      length( unique( sim$population[sim$G,] ))
    ),
    3,
    line = 1
  )

  # Frequencies
  haplotype.frequencies = compute.haplotype.frequencies( sim )
  variant.frequencies = compute.variant.frequencies( sim )

  palette = rainbow( sim$H )
  blank.plot( xlim = c( 1, sim$G ), ylim = c( 0, 1 ) )
  for( i in 1:sim$H ) {
    points( 1:sim$G, haplotype.frequencies[,i], col = palette[i], type = 'l' )
  }
  axis( 1 )
  axis( 2 )
  mtext( "Haplotype frequencies", 2, line = 2, cex = 0.5 )
  mtext( "Generations", 1, line = 2, cex = 0.5 )
  mtext( "Haplotype frequency evolution", 3, line = 1 )
  
  grid()

  palette = rainbow( sim$L )
  blank.plot( xlim = c( 1, sim$G ), ylim = c( 0, 1 ), ylab = "Variant frequencies")
  for( i in 1:sim$L ) {
    points( 1:sim$G, variant.frequencies[,i], col = palette[i], type = 'l' )
  }
  axis( 2 )
  axis( 1 )
  mtext( "Variant frequencies", 2, line = 2, cex = 0.5 )
  mtext( "Generations", 1, line = 2, cex = 0.5 )
  mtext( "Variant frequency evolution", 3, line = 1 )
  grid()

  # LD matrices
  breaks = c( -0.01, seq( from = 0.1, to = 1, by = 0.1 ))
  ld.colours = c( 'grey95', rev( heat.colors( length(breaks)-2 )))

  starting.ld = compute.ld( sim$haplotypes[sim$population[1,],] )^2
  ending.ld = compute.ld( sim$haplotypes[sim$population[sim$G,],] )^2
  image( t( starting.ld ), breaks = breaks, col = ld.colours, x = 1:sim$L, y = 1:sim$L )
  mtext( "LD (r2) at starting generation", 3, line = 1 )
  image( t( ending.ld ), breaks = breaks, col = ld.colours, x = 1:sim$L, y = 1:sim$L )
  mtext( sprintf( "LD (r2) at generation %d", sim$G ), 3, line = 1 )
  
  legend(
    "topleft",
	col = ld.colours[c(2,4,6,8,10)],
	legend = sprintf( "<=%.2f", breaks[c(3,5,7,9,11)]),
	pch = 15,
	bty = 'n'
  )
  
}
