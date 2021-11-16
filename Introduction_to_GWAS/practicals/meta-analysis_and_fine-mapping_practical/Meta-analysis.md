### Meta-analysing the two studies.

#### Conducting fixed-effect meta-analysis
The simplest way to perform meta-analysis is as follows.  We assume:

* that the *true effect* is the same in both studies.
* that the *observed effect* in each study is equal to the *true effect* plus noise.  (The noise is given by the study standard error, and is assumed to be gaussian in the usual way.)

The full name of this is *inverse variance weighted fixed-effect meta-analysis*.  

Here's how it works: form a weighted average of the effect estimates, weighted by the inverse of the variances:

<img src="https://render.githubusercontent.com/render/math?math=b = \sum_i \frac{\beta_i}{\text{se}_i^2}">

Compute also the sum of weights:

<img src="https://render.githubusercontent.com/render/math?math=w = \sum_i \frac{1}{\text{se}_i^2}">

Then the *meta-analysis estimate* is

<img src="https://render.githubusercontent.com/render/math?math=\beta_{\text{meta}} = b * w">

and the *meta-analysis standard error* is

<img src="https://render.githubusercontent.com/render/math?math=\beta_{\text{meta}} = \frac{1}{\sqrt{w}}">

**Note.** There are two good reasons to compute the estimate this way. Firstly, [normal times
normal is
normal](../../Statistical_modelling/Introduction/notes/Normal%20times%20normal%20is%20normal.pdf) -
and if you figure it out, you'll see the above calculation is the same as that lemma.

Secondly, it makes sense: the meta-analysis estimate is a weighted average of the per-study
estimates, and they are weighted by the variance: studies with lots of uncertainty (large variance)
get weighted down while studies with little uncertainty (small variance) get higher weight.

### Running the meta-analsis

Here is a function to run the meta-analysis from our two studies:
```
meta.analyse <- function( beta1, se1, beta2, se2 ) {
  inverse.variances = c( 1/se1^2, 1/se2^2 )
  b = c( beta1, beta2 ) * inverse.variances
  meta.inverse.variance = sum( weights )
  meta.beta = sum( b ) / meta.inverse.variance
  meta.se = sqrt( 1 / meta.inverse.variance )
  
  # This formulation computes a two-tailed P-value
  # as discussed in lectures
  meta.P = pnorm( -abs( meta.beta ), sd = meta.se ) * 2
  return( list(
    meta.beta = meta.beta,
    meta.se = meta.se,
    meta.P = meta.P
  ))
}
```

Try it on some fake data, for example:
```
> meta.analyse( beta1 = 1, se1 = 0.1, beta2 = 1, se2 = 0.1 )
```

Try a few different input values.  You should see:

* The *meta-analysis estimate* should be somewhere between the two input estimates.
* The *meta-analysis standard error* should be smaller than either of the two study standard errors.
* The *meta-analysis P-value* depends on how well the two estimates 'stack up' on each other.  If they are both in the same direction and of roughly comparable size, the P-value will be lower than either of the study P-values.

Here is a way to run the above function across all the data in our study.

**Note.** We will use the [`purrr`](https://purrr.tidyverse.org) function `map_dfr()` to run the above function across all rows and return a data frame.  If you don't have `purrr`, you could use a base R function like `lapply()` instead.

```
meta_analysis_results = map_dfr(
  1:nrow(study1),
  function( i ) {
    c(
      rsid = study1$rsid[i],
      chromosome = study1$chromosome[i],
      position = study1$position[i],
      allele1 = study1$allele1[i],
      allele2 = study1$allele2[i],
      study1_beta = study1$beta[i],
      study1_se = study1$se[i],
      study2_beta = study2$beta[i],
      study2_se = study2$se[i],
      meta.analyse(
        study1$beta[i], study1$se[i],
        study2$beta[i], study2$se[i]
      )
    )
  }
)
```

### Making a forest plot

To make best sense of this the tool for the job is a [*forest plot*](https://en.wikipedia.org/wiki/Forest_plot).  We plot the study estimates and the meta-analysis estimate, along with their confidence/credible intervals, on seperate lines.  Let's do that now.  First here's our utility function to make a blank canvas:

```
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
```

Now let's make a forest plot function:
```

forest.plot <- function(
  betas,
  ses,
  names = c( "Study 1", "Study 2", "Meta-analysis" )
) {
  # y axis locations for the lines.  We separate the meta-analysis line slightly.
  # Note: we assume three betas here: study1, study2, meta-analysis!
  y = c( 3, 2, 0.5 )
  
  # learn a good x axis range by going out 3 standard errors from each estimate:
  xlim = c( min( betas - 3 * ses ), max( betas + 3 * ses ))
  # Also let's make sure to include zero in our range
  xlim[1] = min( xlim[1], 0 )
  xlim[2] = max( xlim[2], 0 )

  # expand the range slightly
  xcentre = mean(xlim)
  xlim[1] = xcentre + (xlim[1] - xcentre) * 1.1
  xlim[2] = xcentre + (xlim[2] - xcentre) * 1.1

  # Give ourselves a big left margin for the row labels
  par( mar = c( 4.1, 6.1, 2.1, 2.1 ))
  blank.plot(
    xlim = xlim,
    ylim = c( range(y) + c( -0.5, 0.5 )),
    xlab = "Estimate and 95% CI"
  )
  
  # Draw the intervals first so they don't draw over the points
  segments(
    x0 = betas - 1.96 * ses, x1 = betas + 1.96 * ses,
    y0 = y, y1 = y,
    col = 'grey'
  )

  # Now plot the estimates
  points(
    x = betas,
    y = y,
    col = 'black',
    pch = 19
  )

  # ... and add labels.  We put them 10% further left than the leftmost point
  # and we right-align them
  text.x = xcentre + (xlim[1] - xcentre) * 1.1
  text(
    x = text.x,
    y = y,
    labels = names,
    adj = 1, # right-adjust
    xpd = NA # this means "Allow text outside the plot area"
  )
  
  # Add an x axis
  axis( 1 )
  
  # add grid lines
  grid()
  
  # add a solid line at 0
  abline( v = 0, col = rgb( 0, 0, 0, 0.2 ), lwd = 2 )
}
```

For example, we could now make a forest plot for the first SNP:
```
  meta = meta.analyse( study1$beta[1], study1$se[1], study2$beta[1], study2$se[1] )
  betas = c( study1$beta[1], study2$beta[1], meta$meta.beta )
  ses = c( study1$se[1], study2$se[1], meta$meta.se )
  forest.plot( betas, ses )
```

<img src="solutions/forest_plot_1st_SNP.png">

**Question.** What does the forest plot look like for the 'top' SNPs, i.e. for those with the lowest P-values.  (Or Bayes factors)?  Plot a few of them and look at them.

**Note.** Forest plots are deceptively simple and you can spend a lot of time tweaking them.  Can you add text to the plot listing the numerical values (estimate, 95% confidence interval, and P-vale?)  E.g. this code:

```
pvalues = pnorm( abs( betas ), sd = ses ) * 2
labels = sprintf(
  "%.2f (%.2f - %.2f) P = %.3g",
  betas,
  betas - 1.96 * ses,
  betas + 1.96 * ses,
  pvalues
)
```
produces the right kind of text.  Then you can add it to the plot in the right place.

