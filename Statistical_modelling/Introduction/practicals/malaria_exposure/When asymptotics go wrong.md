# When asymptotics go wrong

If you've been following the practical so far, you'll have a plot of the parameter estimates from a binomial likelihood for both rows of our table.

Look at the top row of the plot.  This picture is typical for binomial likelhoods when the frequency is low or there's not much data - the asymptotic approximation does not work well.  The likelihood becomes skewed, and the normal approximation starts to fail (here overstating the amount of uncertainty on the left, and understating it on the right.)

To show this more clearly - because this is just a binomial, we can also do an exact computation of the estimate and its confidence interval (this function also computes a P-value against the null model that the frequency is zero) using the `binom.test()` function, like this:
```R
> binom.test( n = 1966, x = 1, p = 0 )

	Exact binomial test

data:  1 and 1966
number of successes = 1, number of trials = 1966, p-value < 2.2e-16
alternative hypothesis: true probability of success is not equal to 0
95 percent confidence interval:
 1.287774e-05 2.830707e-03
sample estimates:
probability of success 
           0.000508647 

```

The estimate is `0.000508647` (which is 1/1966) and the 95% confidence interval is reported as `1.287774e-05` to `2.830707e-03`.  On the other hand the asymptotic method above computes the standard error as `0.0005093123` (see the function output) so we can compute asymptotic estimated confidence intervals as:
```R
theta = 1/1966
se = 0.0005093123
sprintf( "95%% CI = %.4f (%.4f - %.4f)", theta, theta - 1.96 * se, theta + 1.96 * se )
```

This reports `95% CI = 0.0005 (-0.0005 - 0.0015)`.  As the plot indicates, this asymptotic confidence interval overstates uncertainty on the left but understates it on the right.  Similarly, an asymptotic computation of a P-value (for example the Wald test P-value, which just looks at how far the estimate is in the tails of the approximating normal distribution) gets it somewhat wrong - here is a **one-sided test** against the null that ϑ₁ = 0:
```
# compute Wald test p-value using cumulative distribution function of Gaussian:
pnorm( theta, mean = 0, sd = 0.0005093123, lower.tail = F ) 
```

This reports the P-value as `0.15`.  But, since we've observed one C allele, it is *completely impossible* that the frequency is zero!  The true p-value is actually zero, because ϑ₁ = 0 *never* generates any C alleles at all.

On the other hand, when there's lots of data and the frequency is not close to zero (e.g. the second row of the plot) things are fine.  The likelihood becomes a fairly well-rounded citizen (literally so - it becomes approximately normal, i.e. its log has approximately constant curvature).
