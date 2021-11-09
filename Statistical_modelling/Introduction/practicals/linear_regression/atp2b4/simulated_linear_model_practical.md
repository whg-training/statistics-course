## Linear regression for simulated data

### simulating some data
One of the simplest ways to explore statistical models is by simulating.  Let's simulate some data from a simple linear regression and then fit it.

First let's make a gaussian predictor variable:
```
  N = 100
  X = rnorm( N )
```

Now let's make a gaussian outcome variable that is correlated with X. To do this, we'll simply take some multiple of X and add some noise:
```
Y = 0.5*X + rnorm( N, sd = sqrt(0.75) )
```

(I chose that standard deviation because of way variance scales.  If *X* has variance 1, then 0.5 × *X* has variance 0.25.  Therefore we need to add a variable with variance 0.75 to get back to variance 1.  And of course the standard deviation is the square root of the variance.)

Let's see what that looks like

```
library( tidyverse )
data = tibble( X = X, Y = Y )
plot( data$X, data$Y, xlab = "X", ylab = "Y", pch = 19 )
```

**Note.** I am using the [tidyverse](https://www.tidyverse.org) for these examples. If you can't install it (or don't
want to), you can use `data.frame()` instead of `tibble()` above.

### What is linear regression?

See the ![linear regression diagram](https://github.com/whg-gms/statistics-course/raw/main/Statistical_modelling/Introduction/notes/Linear%20regression.pdf)

### fitting the regression

The function `lm()` in R fits a linear regression model.  The syntax is:

```
lm( [outcome] ~ [predictor], data = [data frame] )
```

So let's fit our data:
```
fit = lm( Y ~ X, data = data )
```

**Note.** This fits the linear regression likelihood, finding its maximum.  It does not directly support including a prior distribution.  We'll work with that for now.

If you run this you'll see something like: 

```
> lm( Y ~ X, data = data )

Call:
lm(formula = Y ~ X, data = data)

Coefficients:
(Intercept)            X  
   -0.09646      0.54167  
```

A better summary is provided by (guess what?) `summary()`:

```
fit = lm( Y ~ X, data = data `
summary( fit )
```

In particular it's useful to grab the **coefficients** (i.e. the maximum likelihood parameter values):

```
coeffs = summary( fit )$coeff
print( coeffs[,1:2])
```

**Question**. What do these coefficients mean?  

Let's plot them:

```
# plot the data:
plot(
  data$X,
  data$Y,
  pch = 19,
  xlab = "X",
  ylab = "Y"
)

abline( coef = coeffs[,1], col = 'red', lwd = 2 )
```

That's just the best-fitting line... let's also plot 95% upper and lower estimates of the slope:
```
abline( coef = coeffs[,1] + c( 0, coeffs[2,'Std. Error'] * 1.96 ), col = 'red', lty = 2 )
abline( coef = coeffs[,1] - c( 0, coeffs[2,'Std. Error'] * 1.96 ), col = 'red', lty = 2 )
```

**Note.** The value 1.96 used above comes from the following.  If you work out the area of 95% mass under a standard Gaussian distribution:
```
qnorm( c( 0.025, 0.975 ))
```
you will see it stretches from about -1.96 to 1.96.  That is, a variable that is Gaussian distributed has 95% of its mass within +/- 1.96 standard deviations of the mean.

In the simulation above, the true effect size was 0.5 and the intercept (i.e. the mean of Y) was zero.  Let's plot that as well:
```
abline( coef = c(0,0.5), col = 'grey', lty = 2 )
```
if your simulation is like mine... the fit is 'not bad'.

### What was the likelihood?

The logarithm of the likelihood can be obtained with the `logLik()` function:
```
logLik(fit)
```
**Note.** Although you can exponentiate to get the (non-logged) likelihood, it is very very small. Working in log space is typically better.

**Question.** suppose you go back and simulate again with twice as many data points.  How would the (log)likelihood change?  Would it get bigger or smaller?  Why?

### Regression diagnostics

The linear regression model says that the points should be distributed around the line with the same error distribution - a Gaussian with mean 0 and fixed variance.  Let's make a plot to test that.

First, the fitted residual standard deviation can be obtained with the `sigma` function:
```
residual.sd = sigma(fit)
```

Now let's generate N quantiles of a Gaussian with that standard deviation - these are what we *expect* the data to look like if the model is good:

```
# I use 1:N/(N+1) to get N equally spaced quantiles
expected.residuals = qnorm( 1:N/(N+1), sd = residual.sd )
# You could also use sort( rnorm( N, sd = residual.ld )) here
```

Now let's plot a **qq-plot** comparing the expected to the observed residuals:
```
# two ways to get residuals - calculate directly:
observed.residuals = sort( data$Y - predict(fit) )
# or just get them from the fit object:
observed.residuals = sort( fit$residuals )

plot(
  expected.residuals,
  observed.residuals,
  pch = 19,
  xlab = "Expected residuals",
  ylab = "Observed residuals"
)
abline( a = 0, b = 1, col = 'red' )
```

They look... pretty close to the line?

**note.** Incidentally it is useful to have a function to do this.  How about:
```
lm.qqplot = function( fit ) {
  N = nrow( fit$model )
  residual.sd = sigma(fit)
  expected.residuals = qnorm( 1:N/(N+1), sd = residual.sd )
  observed.residuals = sort( fit$residuals )
  plot(
    expected.residuals,
    observed.residuals,
    pch = 19,
    xlab = "Expected residuals",
    ylab = "Observed residuals"
  )
  abline( a = 0, b = 1, col = 'red' )
}
lm.qqplot(fit)
```

## Fitting to real data

Now go to the [ATP2B4 expression practical](atp2b4_practical.md).
