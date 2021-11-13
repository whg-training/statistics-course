[Up to the table of contents](Introduction.md) - [Back to the previous page](relatedness_pruning.md) - [Forward to the next page](global_analysis.md).

### Computing principal components

If you've reached this page, you should now have downloaded the data, [computed an LD-pruned set of
SNPs](ld_pruning.md), and also [computed a set of largely unrelated samples to work
with](relatedness_pruning.md). Congratulations - you're now ready to compute principal components!

Use plink to compute the PCs:

```
$ plink --vcf chr19-clean.vcf.gz --extract chr19-clean.prune.in --remove related_samples.txt --pca var-wts --out chr19-clean
```

Leave this to run for a minute or so.

**Question**. It's always worth inspecting screen output to check things look right.  For example, has plink excluded the right number of samples we told it to?

Towards the end of the output plink will tell us where it has saved the results. These should be in
the files `chr19-clean.eigenvec` (which stores the actual PCs), `chr19-clean.eigenvec.var` (which
stores the SNP weights or loadings, reflecting how much each SNP contributes to each PC), and
`chr19-clean.eigenval` (which says how much of the overall genotypic variance each PC explains).

#### An aside on the maths.

The maths of PCA (in the way usually applied) works like this. Suppose `X` is a big matrix of
genotypes at the `L` genetic variants (rows) and `N` samples (in columns).

It is typical to always standardise the genotypes at by dividing by each variant (row) in *X* by its standard deviation, i.e. by <img src="https://render.githubusercontent.com/render/math?math=\sqrt{f(1-f)}"> where *f* is the
allele frequency of the variant - and then subtracting the mean.  So do that first.

Now there are two possible square matrices you can make out of *X*:

1. You can form the *relatedness matrix*:

<img src="https://render.githubusercontent.com/render/math?math=R=\frac{1}{L} X^t X">.
  
  *R* is an *N &times; N* matrix i.e. has one row and one column for each
  sample. The motivation is that each entry *r<sub>i,j</sub>* captures the degree of allele sharing
  (covariance) between individual i and j. Also, because of the frequency scaling, sharing of rarer
  variants has greater weight *r<sub>i,j</sub>*.

2. Or you can form the *LD matrix*:

<img src="https://render.githubusercontent.com/render/math?math=Z=\frac{1}{N} X X^t">.

Again because of the scaling to unit variance, entry *z<sub>i,j</sub>* is the LD (correlation) between genotypes at variants i and j, and the diagonal entries are *1*.

If you followed the tutorial this morning you will realise that *relatedness* and *LD* are in some
sense dual to each other. This duality appears in principal components analysis too. Namely:

* The right eigenvectors of *Z* are the *principal component loading vectors* of *X*.

* The right eigenvectors of *R* are (proportional to) the *principal component vectors* (or just
  "principal components") of *X*.

* The two are closely related, since the projection of each sample's genotypes onto the *i*th principal component loading vector, is the *i*th principal component for that sample.

The key property that makes these vectors useful is this: the first loading vector picks out the
direction in genotype space (a linear combinations of SNPs) so that the first principal component
has the maximum possible variance. The second loading vector then picks out the direction in
genotype space *orthogonal to the first* that makes the second principal component have the maximum
possible variance - and so on. These vectors thus pick out *the directions of greatest variance in
the genotype data*.

**Note.** There are two sources of LD that could turn up here. The first derives from genetic drift
or directional selection, which (as we discussed this morning) causes LD between local SNPs on a
chromosome. The second is LD deriving from population structure. This type of LD can occur between
SNPs anywhere in the genome - for example, it would be seen between any SNPs that differ in
frequency between sub-populations. Because we have ld-thinned our data, it is primarily this type
of population structure-driven LD that will be picked up by our principal components analysis.

**Computation in practice**. Typically *Z* is huge, while *R* is somewhat smaller.  So tools like
`plink` compute *R*, use this to compute the principal components, and then compute loadings using
a second pass through the data.

You can easily do this in R as well, for example:
```
X = (load data here, convert to a matrix and standardise...)
L = nrow(X)
R = (1/L) * t(X) %*% X
PCs = eigen(R)$vectors
plot( PCs[,1], PCs[,2], xlab = "PC 1", ylab = "PC 2" )
# etc.
```

The `eigen()` step will take quite a while on any real dataset, but is generally manageable up to around a few thousand samples.

**Computing in big cohorts**. While the above works in many studies, very large cohorts such as the
[UK Biobank[(https://www.nature.com/articles/s41586-018-0579-z) typically require other methods such as [flashPCA](https://github.com/gabraham/flashpca) that avoid computing either of the matrices directly.

**Note.** If all this hasn't melted your brain, [try reading
Gil McVean's paper on genealogies and PCA](https://doi.org/10.1371/journal.pgen.1000686).

### Plotting the principal components

Now let's load the PCs and plot them:

In R/RStudio:

```
pcs = read.table( "chr19-clean.eigenvec" )
View(pcs)
plot( pcs[,3], pcs[,4], pch = 19, xlab = "PC1", ylab = "PC2" )
```

**Question**. What do you see?  Is there any obvious structure in the dataset?  

**Note.** I used `pch = 19` above because that gives us solid dots instead of open circles.  See [here for a diagram of available shapes in R](https://r-graphics.org/recipe-scatter-shapes).

We might also want to plot more than just the top two PCs. Let's plot all pairs of the top 5 PCs:

In R/RStudio:

```
colnames(pcs)[3:7] = c( "PC1","PC2","PC3","PC4","PC5" )
pairs( pcs[, 3:7], pch = 19 )
```

**Hint**: In Rstudio, you can click the 'zoom' button above above the plot to get a bigger version.

Our dataset contains samples with different ethnicities that in our dataset are recorded in the
first column of the file (plink calls this a 'family ID' but here we've used it to record sampled
ethnicities). Make a list of them using `table()`:

```
table( pcs[,1] )
```

Different ethnic groups might be genetically distinct, so they might separate on the PCs. Let's
replot colouring points by their ethnicity. We'll do this in three steps:

1. Turn ethnicities into integers that we can pass to the col argument of `plot()`.
2. Plot using the colours
3. Add a legend to the plot so we can see what is what.

**Note.** The `lower.panel` and `xpd` arguments below are used to tweak the appearance of the plot and how the legend is plotted.  See `?pairs` and `?par` if you want more information, but don't worry about these for now.

In R/RStudio:

```
colnames(pcs)[1] = "group"
palette = c( CAN = "red2", FAN = "green2", JAN = "blue2", WAN = "yellow3" )
pcs$colour = palette[ pcs$group ]
pairs( pcs[,3:7], col = pcs$colour, pch=20, lower.panel = NULL )
legend(
  0.1, 0.5,
  col = palette,
  legend = names( palette ),
  pch=20,
  xpd=NA   # this is used to allow drawing outside the plot.
)
```

**Question**. Do the principal components reflect ethnicities in this dataset?  Which ethnicities are separated by which PCs? 

**Question**. Do any samples look especially genetically distinct from the other samples? What might be the reason for this?  Is the same set of samples outlying on all the PCs?  Use your R skills to identify these samples.  

### Plotting SNP weights/loadings

But we are not done yet! If we want to look at population structure, we generally want PCs that
represent genome-wide variation in allele frequencies (as opposed to variants clustered in specific
regions.) To check this we should also plot the loadings.

To do this, examine the `chr19-clean.eigenvec.var` file that plink produced. We'll load it
into R and plot loadings across the chromosome.

In R/RStudio:

```
loadings = read.table("chr19-clean.eigenvec.var")
View(loadings)
```

Columns 3-22 of this file represent the loadings on PCs 1-20, respectively.  Let's plot loadings for the first 5 PCs.  Here we can use the `layout()` function to make a plot with multiple rows.  The `par( mar = ...)` command adjusts the plot margins.)

In R/RStudio:

```
layout( matrix( 1:5, ncol = 1 ))
par( mar = c( 1, 4, 1, 2 ))
for( i in 1:5 ) {
  plot(
    1:nrow(loadings),

    # we plot the absolute value of the loading, so they all go up
    abs( loadings[,i+2] ),  

    ylab = paste("PC ", i, sep = "" ),
    ylim = c( 1, 10 )
  )
}
```

**Question**. Do any SNPs stand out as having high loadings for particular PCs?  Which PCs?  Look back at your plot of principal components.  Do these PCs pick out particular clusters of samples?

**Question**. Go back and compute PCs without excluding related samples.  Then re-load and plot the PCs and loadings.  What drives the top PCs now?  What do the loadings look like?

### Plotting samples against a global dataset

A last step in many analyses is to check the ancestry of samples against a global dataset - to check for example that no unexpected ancestries are included.  We'll do that [on the next page](global_analysis.md).
