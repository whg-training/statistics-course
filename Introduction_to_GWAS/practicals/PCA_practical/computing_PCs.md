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

If you want to know more about what is being computed and how - see [this aside on the maths of principal components analysis](the_maths.md).

And if all this isn't enough [try reading
this seminal paper on sample genealogies and PCA](https://doi.org/10.1371/journal.pgen.1000686).


### Plotting the principal components

Now let's load the PCs and plot them:

In R/RStudio:

```
pcs = read.table( "chr19-clean.eigenvec" )
colnames(pcs)[1:7] = c( "group", "ID", "PC1","PC2","PC3","PC4","PC5" )
View(pcs)
plot( pcs$PC1, pcs$PC2, pch = 19, xlab = "PC1", ylab = "PC2" )
```

**Question**. What do you see?  Is there any obvious structure in the dataset?  

**Note.** I used `pch = 19` above because that gives us solid dots instead of open circles.  See [here for a diagram of available shapes in R](https://r-graphics.org/recipe-scatter-shapes).

We might also want to plot more than just the top two PCs. Let's plot all pairs of the top 5 PCs:

In R/RStudio:

```
pairs( pcs[, 3:7], pch = 19 )
```

**Hint**: In Rstudio, you can click the 'zoom' button above above the plot to get a bigger version.

Our dataset contains samples with different ethnicities that in our dataset are recorded in the
first column of the file (plink calls this a 'family ID' but here we've used it to record sampled
ethnicities). Make a list of them using `table()`:

```
table( pcs$group )
```

Different ethnic groups might be genetically distinct, so they might separate on the PCs. Let's
replot colouring points by their ethnicity. We'll do this in three steps:

1. Turn ethnicities into integers that we can pass to the col argument of `plot()`.
2. Plot using the colours
3. Add a legend to the plot so we can see what is what.

**Note.** The `lower.panel` and `xpd` arguments below are used to tweak the appearance of the plot and how the legend is plotted.  See `?pairs` and `?par` if you want more information, but don't worry about these for now.

In R/RStudio:

```
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
colnames(loadings)[1:4] = c( "chromosome", "rsid", "alleleA", "alleleB" )
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
    abs( loadings[,i+4] ),  
    ylab = paste("PC ", i, sep = "" ),
    ylim = c( 1, 10 ),
    pch = 19
  )
  legend( "topleft", legend = sprintf( "Loadings for PC %d", i ), bty = 'n' )
}
```

**Question**. Do any SNPs stand out as having high loadings for particular PCs?  Which PCs?  Look back at your plot of principal components.  Do these PCs pick out particular clusters of samples?

**Question**. Go back and compute PCs without excluding related samples.  E.g. something like this, without the `--remove` option:

```
$ plink --vcf chr19-clean.vcf.gz --extract chr19-clean.prune.in --pca var-wts --out chr19-all_samples
```

Then re-load and plot both the PCs and loadings.  What drives the top PCs now?  What do the loadings look like?

### Plotting samples against a global dataset

A last step in many analyses is to check the ancestry of samples against a global dataset - to check for example that no unexpected ancestries are included.  We'll do that [on the next page](global_analysis.md).
