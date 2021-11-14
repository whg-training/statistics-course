[Up to the table of contents](Introduction.md) - [Back to the page on LD pruning](ld_pruning.md) - [Forward to start computing PCs](computing_PCs.md)

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

#### Picking samples to exclude

Before we compute PCs, we'll pick one sample from each closely-related pair and exclude it. For
now we'll just do this simply by picking the second sample; a more refined approach might look at
genotyping performance across the genome and exclude the sample with greater missingness. (Better
yet, you could try to pick the fewest number of samples that remove all close relationships, as
this can sometimes substantially improve the resulting sample size). We'll store the results in a
file that we can tell plink about in later steps.

In RStudio:

```
# put a line on our plot so we can see what we'll exclude
abline( v = 0.2, col = 'red' )

# Write a list of excluded samples
exclusions = ibd[ ibd$PI_HAT > 0.2, c('FID2','IID2') ]
write.table( exclusions, file="related_samples.txt", col.names = F, row.names = F, quote = F )
```

Q. How many samples will be excluded?

**Note**: the same sample might appear twice in the exclusions data frame.  This would happen if they are closely related to more than one other sample.  Are there any samples like that?  You can use `unique(exclusions)` to get a list of unique samples to exclude.

#### Computing principal components

When you're ready, go [here](computing_PCs.md) to start the principal component analysis proper.

