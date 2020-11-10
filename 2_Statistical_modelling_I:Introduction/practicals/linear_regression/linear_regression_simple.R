X = read.table( "linear_regression_simple.tsv", hea=T, as.is=T )
l = lm( outcome ~ A + B, data = X )

# Look at coefficients:
summary(l)$coeff

# Plot outcome versus predictor
plot( X$A, X$outcome, pch = 19, bty = 'n' )
abline( l, col = 'red' )
grid()

sigma = sigma(l)  # get residual standard deviation
# get equally spaced quantiles of the normal distribution
expected = qnorm(
  (1:nrow(X))/(nrow(X)+1),  # equally spaced points in the unit interval,
  sd = sigma
)

plot(
  expected,                #actual normal quantiles
  sort( residuals(l) ),    #observed residuals sorted
  pch = 19
)
abline( a = 0, b = 1, col = 'red' ) # diagonal line
grid()


# Try a model with non-gaussian residuals
l = lm( exp(outcome) ~ A + B, data = X )
plot( X$A, exp(X$outcome), pch = 19, bty = 'n' )
abline(l, col = 'red' )

sigma = sigma(l)  # get residual standard deviation

# get equally spaced quantiles of the normal distribution
expected = qnorm(
  (1:nrow(X))/(nrow(X)+1),  # equally spaced points in the unit interval,
  sd = sigma
)

plot(
  expected,                #actual normal quantiles
  sort( residuals(l) ),    #observed residuals sorted
  pch = 19
)
abline( a = 0, b = 1, col = 'red' ) # diagonal line
grid()

