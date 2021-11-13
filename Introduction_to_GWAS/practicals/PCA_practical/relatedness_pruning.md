[Up to the table of contents](Introduction.md) - [Back to the previous page](ld_pruning.md) - [Forward to the next page](computing_PCs.md)

### IBD pruning of samples and identification of close relationships

We want our top PCs to reflect the relatedness structure across the majority of samples in our
GWAS dataset. Although all samples in our dataset are nominally unrelated, a few duplicated or
related samples may have slipped in through the sampling process or through sample handling. As
described in the lectures, this can seriously affect principal components. We'll therefore first
identify and remove any close relationships before computing PCs.

Let's use plink to compute the relatedness between samples. A simple way to do this would be to use
the relatedness matrix that we will construct to compute PCA. However, here we'll use genome-wide
estimates of identity by descent (IBD) instead. To do this we use the following command:

```
$ plink --vcf chr19-clean.vcf.gz --genome gz --out chr19-clean --extract chr19-clean.prune.in
```

Here, the `--genome` option tells plink to compute genome-wide IBD estimates (as this file is rather large, we ask plink to store it in a gzipped file).  We've used the `--extract` option to tell plink to only use the LD-pruned set of SNPs we computed above.  Further details of the options for computing relatedness and file formats are described [on the plink documentation page](https://www.cog-genomics.org/plink2/ibd).

The output file from the command above is in `chr19-clean.genome.gz`. Let's load that into `R`
and have a look at it.

In R/RStudio:

```
ibd <- read.table("chr19-clean.genome.gz", hea=T, as.is=T)
View(ibd)
hist( ibd$PI_HAT, breaks = 100 )
```

You can also zoom in along the y-axis to see any close relationships:

```
hist( ibd$PI_HAT, breaks = 100, ylim = c(0,1000) )
```

**Question**. What is `PI_HAT`? Can you figure it out from [the documentation
page](https://www.cog-genomics.org/plink2/ibd)?

**Question**. Which samples have the closest relationship in the dataset?  What is the average relationship between samples in this dataset?

**Question**. In the histogram, there's a tail or 'bump' of relationships up to about 0.2.  What could this represent?

**Note.** The above is one way to compute relatedness estimates - many other methods are available.
For example, [KING](https://www.kingrelatedness.com) is a popular choice.  More computationally demanding methods also exist that can [identify actual segments inherited IBD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7553009/).

#### Computing principal components

When you're ready, go [here](computing_PCs.md) to start the principal component analysis proper.

