## The frequency of O blood group

In this section we will estimate the frequency of O blood group (as encoded by the [common
loss-of-function deletion in the ABO
gene](https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs8176719), in multiple populations.

#### An aside on the ABO blood group system.

The ABO blood group system was the first genetic polymorphism discovered in humans - [long before
the structure of DNA was solved](https://www.ncbi.nlm.nih.gov/books/NBK2267/). It was discovered by
studying agglutination patterns of red cells in serum from other individuals and is of course of
extreme relevance to blood transfusion. Beyond this, however, ABO is an interesting gene; the A/B
split seems to have been preserved [under balancing
selection](https://www.pnas.org/content/109/45/18493) across primates, while the O mutation itself
may be a recurrent polymorphism. ABO is also [one of the most pleiotropic loci in the human
genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5207801/) and there is evidence it is
[associated with risk of SARS-CoV-2
infection](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7594382/). 

But I was interested because it is also [protective against
malaria](https://www.nature.com/articles/s41467-019-13480-z). (The data here comes from that paper).

### Loading the data
First, let's load the data.

**Note.** I am using the [tidyverse](https://www.tidyverse.org) for these examples. If you don't
have it installed, either install it or use base R (e.g. `read.csv()`) instead.

```
library( tidyverse )
data = read_csv( "o_bld_grp.csv" )
```

**Question**.  What countries are in the data?  What ethnic groups?

The O blood group data is in the column called `O.bld.grp`. A `1` in this column means the
individual has O blood group (which happens if they have two copies of the above rs8176719
deletion). A `0` generally means they will have either A, AB, or B blood group depending on the
alleles they carry. (There are actually a few other mutations that cause loss of function of `ABO`,
but we're ignoring them here.)

**Question**. Table this data for each population or ethnicity. Can you estimate the frequency in
each population? In each ethnic group?


### Building a model

You probably estimated the allele frequency as: the number of O individuals divided by the total
number of individuals.  If so, good work!

**Question.** Do you trust these estimates? How much? Do you trust your estimates for all the populations the same
amount? Why?

In this course we are all about handling uncertainty. What we want to do is, not just generate a point estimate of the
frequency, but also quantify the uncertainty we have about it.

To do this let's focus on a single population first - say Tanzania:

```
w = which( data$country == "Tanzania" )
> table( data$country[w], data$O.bld.grp[w] )
          
            0  1
  Tanzania 61 41
```

Counts like this can be modelled well using a [binomial distribution](../../notes/Distributions%20cheatsheet.pdf).

**Question**. Suppose we use a binomial distribution to model these counts.  What assumptions are we making?

The binomial distribution takes two parameters: `n`, the number of 'trials' (i.e. samples), and &theta;, the frequency.
In our setting we imagine &theta; to be the 'true' frequency, which is what we want to infer.

The basic inference formula (Bayes rule) is:

<img src="https://render.githubusercontent.com/render/math?math=P(\theta=x|\text{data}) = \frac{P(\text{data}|\theta=x) \cdot P(\theta=x)}{P(\text{data})}">

For the moment let us **assume the prior term *P(&theta;)* is uniform (i.e. ignore it).  


**Note.** As in the [probability cheatsheet](../../notes/Probability%20cheatsheet.pdf) I should stress that all
probability is conditional.  Foir example, 

###

## A note on binomial assumptions

Above I asked what assumptions we make in choosing a binomial distribution. We are making quite a few:

* We are assuming that the total number of samples (e.g. 102 in the case of Tanzania) was known beforehand. (This could
  be violated by sampling. For example this would be violated if we had chosen to up-sample ethnic groups with higher O
  blood group frequency.)
  
* We are assuming that the data points from different individuals can be treated as independent. (This assumption would
  be violated, for example, if we sampled within families, since they share DNA.)




