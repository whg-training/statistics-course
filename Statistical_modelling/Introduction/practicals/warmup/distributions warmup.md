### Plotting distributions


**Question.** Pick a distribution from the [distributions cheatsheet](../notes/Distributions
cheatsheat.pdf) and explore it using the ["distributions zoo"](https://ben18785.shinyapps.io/distribution-zoo/).

**Question.** Pick a distribution from the [distributions cheatsheet](../notes/Distributions
cheatsheat.pdf) and plot its pdf and its cdf for a few different parameter values using R.

**Challenge question.** Plot a multivariate normal distribution density using using R.

**Hint.** This can be done using the `mvtnorm` package and the `ellipse` package.  Here is one way using base R graphics - experiment with different :

```

plot_bivariate_gaussian <- function(
  mean = c( 0, 0 ),
  covariance = matrix(
    c(
      1.0, 0.5,
      0.5, 1.0
    ),
    nrow = 2
  )
) {
  library( mvtnorm )
  library( ellipse )

  # choose a mean and covariance for our MVN
  sigma = covariance

  # generate some samples
  samples = as_tibble( rmvnorm( 1000, mean = mean, sigma = sigma ))
  colnames(samples) = c( "x", "y" )

  # Generate a blank plot
  range = c( -5, 5 )
  plot( 0, 0, col = 'white', xlim = range, ylim = range, bty = 'n' )
  points(
    samples$x, samples$y,
    col = rgb( 0, 0, 0, 0.2 ), # use transparent black
    pch = 19, # small round filled dots
    cex = 0.75 # make points 75% of default size
  )
  for( quantile in c( 0.95, 0.75, 0.5, 0.25, 0.05 )) {
      e = ellipse( sigma, centre = mean, level = quantile )
      points( e[,1], e[,2], type = 'l', col = 'red', lwd = 2 )
  }
  grid()
  abline( h = 0, col = rgb( 0, 0, 0, 0.2 ) )
  abline( v = 0, col = rgb( 0, 0, 0, 0.2 ) )
}

plot_bivariate_gaussian(
  mean = c( 0, 0 ),
  covariance = matrix( c( 1, 0, 0, 1 ), nrow = 2 )
)

```

What happens as you change the mean and/or the covariance?  (Note: the covariance matrix must be **symmetric**, i.e it should have the same value in the upper-right and lower-left entries.  And the diagonal entries must also be positive: they are the marginal variances of the two variables.

### Coin tosses

I toss a coin 100 times.  

**Question.** Covid-19 lateral flow tests currently in use [are thought to be pretty
accurate](https://www.ox.ac.uk/news/2020-11-11-oxford-university-and-phe-confirm-lateral-flow-tests-
 show-high-specificity-and-are). According to that article, the false positive rate is around 0.32%:

<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{not infected}) = 0.0032">

and the specificity is around 80%:

<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{infected}) = 0.8">

Currently [about 0.4% of people in Oxfordshire are infected](https://phdashboard.oxfordshire.gov.uk):

<img src="https://render.githubusercontent.com/render/math?math=P(\text{infected}) = 0.004">

Suppose you test positive.  How worried are you that you have COVID?
