## Estimating an allele frequency

In this section we will estimate the frequency of O blood group (as encoded by the [common
loss-of-function deletion in the ABO gene,
rs8176719](https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs8176719) in multiple
populations. 


### Loading the data
First, let's load the data.

**Note.** I am using the [tidyverse](https://www.tidyverse.org) for these examples. If you don't
have it installed, either install it or use base R (e.g. `read.csv()`) instead.

```
library( tidyverse )
data = read_csv( "../o_bld_grp/o_bld_grp.csv" )
```

**Question**.  What countries are in the data?  What ethnic groups?

The O blood group data is in the column called `O.bld.grp`. A `1` in this column means the
individual has O blood group (which happens if they have two copies of the rs8176719 deletion). A
`0` generally means they will have either A, AB, or B blood group depending on the alleles they carry.

**Question**. Table this data for each population or ethnicity. Can you estimate the frequency in
each population? In each ethnic group?

### Building a model

You probably estimated the allele frequency as: the number of O individuals divided by the total
number of individuals - a good choice.

But in this course we are all about handling uncertainty, and these point estimates aren't good
enough for this.  Instead we need a statistical model to quantify uncertainty.

**Question.** Focus on a single population first, for example Tanzania. 
```
> w = which( data$country == "Tanzania" )
> table( data$country[w], data$O.bld.grp[w] )
          
             0   1
  Tanzania 438 368

```

What distribution should we use to model this data?







## A note on the ABO blood group system.

The ABO blood group system was the first genetic polymorphism discovered in humans - [long before
the structure of DNA was solved](https://www.ncbi.nlm.nih.gov/books/NBK2267/). It was discovered by
studying agglutination patterns of red cells in serum from other individuals and is of course of
extreme relevance to blood transfusion. Beyond this, however, ABO is an interesting gene; the A/B
split seems to have been preserved [under balancing
selection](https://www.pnas.org/content/109/45/18493) across primates, while the O mutation itself
may be a recurrent polymorphism. ABO is also [one of the most pleiotropic loci in the human
genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5207801/) and there is evidence it is
[associated with risk of SARS-CoV-2
infection](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7594382/). (But I was interested because it
is [protective against malaria](https://www.nature.com/articles/s41467-019-13480-z).)
