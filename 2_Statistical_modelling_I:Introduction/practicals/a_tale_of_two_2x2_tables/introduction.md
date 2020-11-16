# A tale of two 2x2 tables

Welcome to (possibly) the world's first *choose your own statistical adventure*!

In this practical we are going to study the following two two-by-two tables, which come from two
published papers about genetic variants that affect malaria susceptibility:

```R
                       non-O      O                                          rs60822373   rs60822373
                       blood gp.  blood gp.                                  G allele     C allele
            controls:  3420       3233               unexposed populations:  1965         1
severe malaria cases:  3925       2738         malaria-exposed populations:  707          17
```

The left table contains data from <https://doi.org/10.1038/s41467-019-13480-z>, and shows that O blood group is at
lower frequency in severe malaria cases than in the general population, consistent with a protective effect of O blood
group.

The right table comes from <https://science.sciencemag.org/content/348/6235/711>, and shows that the rs60822373 'C'
alleles is at higher frequency in malaria-exposed population (here sub-Saharan Africans) than in non-exposed
populations (European-ancestry individuals). rs60822373 encodes the Cromer blood group, so this is consistent with a
protective effect of the Cromer blood group on malaria. (In this right table, the data actually comes from the [1000
Genomes Project phase I](http://www.nature.com/nature/journal/v491/n7422/full/nature11632.html)).

If you run a statistical test on either table you will see that these are highly statistically significant differences
in frequencies between the two rows of each table, so these differences are not just due to sampling.

But *something is wrong with one of these tables*.  Your mission is to find out what!

*QN*.  Maybe you can already see what is wrong?  Give it some thought before we proceed.

## Modelling 2x2 tables

To quantify the effects in both tables we need a model. In both tables, the data was collected based on the row labels
(i.e. case/control or population of collection) so the most straightforward model is of *binomial sampling in rows*.
In other words, for a table:
```

           X    Y        (total)
unexposed  a    b        (N₁ = a+b)
  exposed  c    d        (N₂ = c+d)
```
We model the counts of `Y` as:
```  
  b ~ binomial( n = N₁, p = ϑ₁ )
  d ~ binomial( n = N₂, p = ϑ₂ )
```
where ϑ₁ and ϑ₂ are frequencies in the two rows.

*Exercise* You can implement this model easily using the code we have already written - see [this file](../loglikelihoods/table.ll.R) for an implementation.

However, we'd really like to parameterise the model so that it contains an *effect size parameter* - a parameter that
reflects how different the two row frequencies are.  We will use the *log odds ratio* to do this.  In the above table, the observed log odds ratio is the quantity:

![\text{log} \text{OR}=\log\left(\frac{a d}{b c}\right)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Ctext%7Blog%7D+%5Ctext%7BOR%7D%3D%5Clog%5Cleft%28%5Cfrac%7Ba+d%7D%7Bb+c%7D%5Cright%29)

There are a couple of reasons to focus on log odds ratio here.  One of them is mathematical convenience.  For a probability p, the odds is defined as p/(1-p).  The odds lies in [0,∞] and thus the log-odds lies in [-∞,∞].  It is generally useful to have parameters living in an unbounded space.

Thus, transforming to log odds allows us to put our parameters on the real line.

Another way to think of this is through the inverse function, which maps log odds back to probability space.  This function is known as the *logistic* function and defined by:

![\text{logistic}(x) = \frac{e^x}{1+e^x}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Ctext%7Blogistic%7D%28x%29+%3D+%5Cfrac%7Be%5Ex%7D%7B1%2Be%5Ex%7D)

Let's plot it:
```R
logistic <- function( x ) {
	exp(x) / ( 1 + exp(x) )
}

x = seq( from = -10, to = 10, by = 0.01 )
plot( x, logistic(x), type = 'l', bty = 'n' )
grid()
```

![logistic function](solutions/logistic.png)

You can see that the logistic function maps the real line (log-odds space, x axis) to the unit interval (y axis, probability space).  It is smooth and tails off to being essentially flat outside around [-10,10].

*Exercise* prove that `logistic()` is the inverse of the log odds (substitute one expression into the other and simplify). Alternatively, prove this to yourself by implementing both function in R.

*Exercise* (advanced) what is the slope of the logistic function at x=0?  (Hint: you need to compute the derivative.)

Now that we are working in log-odds space (on the real line), this frees us up to express our log odds as a *baseline log odds* (applying to both rows of the table) plus an *additional log odds conferred by the second row*.  Let's call these parameters μ and β.  We are thus writing:

![\text{log} \text{odds}(\theta_1)=\mu
\quad\text{and}\quad
\text{log} \text{odds}\left(\theta_2\right) = \mu+\beta
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Ctext%7Blog%7D+%5Ctext%7Bodds%7D%28%5Ctheta_1%29%3D%5Cmu%0A%5Cquad%5Ctext%7Band%7D%5Cquad%0A%5Ctext%7Blog%7D+%5Ctext%7Bodds%7D%5Cleft%28%5Ctheta_2%5Cright%29+%3D+%5Cmu%2B%5Cbeta%0A)

Or equivalently:

![\theta_1 = \text{logistic}(\mu)
\quad\text{and}\quad
\theta_2 = \text{logistic}(\mu+\beta)
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Ctheta_1+%3D+%5Ctext%7Blogistic%7D%28%5Cmu%29%0A%5Cquad%5Ctext%7Band%7D%5Cquad%0A%5Ctheta_2+%3D+%5Ctext%7Blogistic%7D%28%5Cmu%2B%5Cbeta%29%0A)

If you follow this maths through you will see the essential fact: parameter β is equal to the log odds ratio, 

![\beta = \log \left(
\frac{\theta_1}{1-\theta_1}
/
\frac{\theta_2}{1-\theta_2}
\right)
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbeta+%3D+%5Clog+%5Cleft%28%0A%5Cfrac%7B%5Ctheta_1%7D%7B1-%5Ctheta_1%7D%0A%2F%0A%5Cfrac%7B%5Ctheta_2%7D%7B1-%5Ctheta_2%7D%0A%5Cright%29%0A)

The reparameterised model for our table is now:

![\left(
\begin{matrix}
1-\text{logistic}(\mu) & \text{logistic}(\mu) \\
1-\text{logistic}(\mu+\beta) & \text{logistic}(\mu+\beta)
\end{matrix}
\right)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cleft%28%0A%5Cbegin%7Bmatrix%7D%0A1-%5Ctext%7Blogistic%7D%28%5Cmu%29+%26+%5Ctext%7Blogistic%7D%28%5Cmu%29+%5C%5C%0A1-%5Ctext%7Blogistic%7D%28%5Cmu%2B%5Cbeta%29+%26+%5Ctext%7Blogistic%7D%28%5Cmu%2B%5Cbeta%29%0A%5Cend%7Bmatrix%7D%0A%5Cright%29)

or fully expanded:
![\left(
\begin{matrix}
\frac{1}{1+e^\mu}
&
\frac{e^\mu}{1+e^\mu}
\\
\frac{1}{1+e^{\mu+\beta}}
&
\frac{e^{\mu+\beta}}{1+e^{\mu+\beta}}
\end{matrix}
\right)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cleft%28%0A%5Cbegin%7Bmatrix%7D%0A%5Cfrac%7B1%7D%7B1%2Be%5E%5Cmu%7D%0A%26%0A%5Cfrac%7Be%5E%5Cmu%7D%7B1%2Be%5E%5Cmu%7D%0A%5C%5C%0A%5Cfrac%7B1%7D%7B1%2Be%5E%7B%5Cmu%2B%5Cbeta%7D%7D%0A%26%0A%5Cfrac%7Be%5E%7B%5Cmu%2B%5Cbeta%7D%7D%7B1%2Be%5E%7B%5Cmu%2B%5Cbeta%7D%7D%0A%5Cend%7Bmatrix%7D%0A%5Cright%29)

Summarising, the log-odds ratio β is a real number which encodes the difference in frequency between the two rows of the table.
