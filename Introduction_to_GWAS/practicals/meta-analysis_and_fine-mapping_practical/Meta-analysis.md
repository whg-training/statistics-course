### Meta-analysing the two studies.

The simplest way to perform meta-analysis is as follows.  We assume:

* that the *true effect* is the same in both studies.
* that the *observed effect* in each study is equal to the *true effect* plus noise.  (The noise is given by the study standard error, and is assumed to be gaussian in the usual way.)

The full name of this is *inverse variance weighted fixed-effect meta-analysis*.  

Here's how it works: form a weighted average of the effect estimates, weighted by the inverse of the variances:

<img src="https://render.githubusercontent.com/render/math?math=b = \frac{\beta_1}{\text{se}_1^2} \plus \frac{\beta_2}{\text{se}_2^2}">

Compute also the sum of weights:

<img src="https://render.githubusercontent.com/render/math?math=w = \frac{1}{\text{se}_1^2} + \frac{1}{\text{se}_2^2}">

Then the *meta-analysis estimate* is

<img src="https://render.githubusercontent.com/render/math?math=\beta_{\text{meta}} = b/w">

and the *meta-analysis standard error* is

<img src="https://render.githubusercontent.com/render/math?math=\beta_{\text{meta}} = \sqrt(w)">

**Note.** Why is it computed this way?

Firstly, [normal times normal is
normal](../../Statistical_modelling/Introduction/notes/Normal%20times%20normal%20is%20normal.pdf) -
and if you figure it out, you'll see the above calculation is the same as that lemma.

Secondly, it makes sense: the meta-analysis estimate is a weighted average of the per-study estimates, and they are weighted by the variance: studies with lots of uncertainty (large variance) get weighted down while studies with little uncertainty (small variance) get higher weight.
