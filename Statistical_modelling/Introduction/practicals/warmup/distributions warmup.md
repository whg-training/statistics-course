## Probability distribution questions

The following can be solved using the [distributions cheatsheet](../../notes/Distributions%20cheatsheet.pdf) and the relevant functions in R.

### Exploring distributions

**Question.** Pick a distribution from the [distributions cheatsheet](../../notes/Distributions%20cheatsheet.pdf) and explore it using the ["distributions zoo"](https://ben18785.shinyapps.io/distribution-zoo/).

**Question.** Pick a distribution from the [distributions cheatsheet](../../notes/Distributions%20cheatsheet.pdf) and plot its pdf and its cdf for a few different parameter values using R.

###

**Question.** What is the mean and variance of a binomial variable (like a coin toss) with parameters *n* (number of trials) and *p* (success probability?)

**Hint.** A binomial variable is a sum of *n* Bernoulli variables (i.e. single coin tosses).  And the [variance of a sum of independent variables is the sum of the variances](https://en.wikipedia.org/wiki/Variance#Basic_properties) (and similarly for the mean).  So you only need to work it out for a single trial and add up.

### Using distribution functions

**Question.** I toss a coin 100 times.  What's the chance I see at least 60 heads?

### Multivariate distributions.

**Question.** Plot a multivariate normal distribution density using using R.

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

  # generate some samples
  samples = as.data.frame( rmvnorm( 1000, mean = mean, sigma = covariance ))
  colnames(samples) = c( "x", "y" )

  # Generate a blank plot
  range = c( -5, 5 )
  plot( 0, 0, col = 'white', xlim = range, ylim = range, bty = 'n' )

  # Plot the samples
  points(
    samples$x, samples$y,
    col = rgb( 0, 0, 0, 0.2 ),  # use transparent black
    pch = 19,                   # small round filled dots
    cex = 0.75                  # make points 75% of default size
  )

  # Plot quantiles as contours
  for( quantile in c( 0.95, 0.75, 0.5, 0.25, 0.05 )) {
      e = ellipse( covariance, centre = mean, level = quantile )
      points( e[,1], e[,2], type = 'l', col = 'red', lwd = 2 )
  }
  # draw some grid lines etc.
  grid()
  abline( h = 0, col = rgb( 0, 0, 0, 0.2 ) )
  abline( v = 0, col = rgb( 0, 0, 0, 0.2 ) )
}

plot_bivariate_gaussian(
  mean = c( 0, 0 ),
  covariance = matrix( c( 1, 0.5, 0.5, 1 ), nrow = 2 )
)

```

Another way to do this would use `ggplot2`'s `geom_bin2d()` to show a 2d histogram:
```
library( ggplot2 )
# (using samples as in the functino above)
ggplot2( data = samples, mapping = aes( x = x, y = y )) + geom_bin2d()
```

What happens as you change the mean and/or the covariance?  (Note: the covariance matrix must have positive diagonal entries: these are the variances of *x* and *y*.  And it must be **symmetric**, i.e it should have the same value in the upper-right and lower-left entries.  This value is the *covariance* between *x* and *y*.)

What do you notice about the contours?
