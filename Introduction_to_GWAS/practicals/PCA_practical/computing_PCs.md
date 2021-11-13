[Up to the table of contents](Introduction.md)

[Back to the previous page](relatedness_pruning.md)

## Computing principal components

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

### An aside on the maths.

The maths of PCA works like this.  Suppose `X` is a big matrix of genotypes at the `L` genetic variants (rows) and `N` samples (in columns).  Then:

* First, standardise the genotypes at each variant by dividing by <img src="https://render.githubusercontent.com/render/math?math=\sqrt{f(1-f)}">, where *f* is the frequency, and then subtracting the mean.
* Second, compute the big <img src="https://render.githubusercontent.com/render/math?math=N\times N"> matrix <img src="https://render.githubusercontent.com/render/math?math=R = \frac{1}{L} X^t X">.  The motivation is that each entry r<sub>i,j</sub> captures the degree of allele sharing (covariance) between individual i and j, with sharing of rarer variants having greater weight.
* The **principal components** are the entries of the (right) eigenvectors of *R*.


# Plot the principal components

Let's load the PCs and plot them:

In RStudio:

```
pcs = read.table( "chr19-clean.eigenvec" )
View(pcs)
plot( pcs[,3], pcs[,4], xlab = "PC1", ylab = "PC2" )
```

**Question**. What do you see?  Is there any obvious structure in the dataset?  

We might also want to plot more than just the top two PCs.  Let's plot all pairs of the top 5 PCs:

In RStudio:

```
colnames(pcs)[3:7] = c( "PC1","PC2","PC3","PC4","PC5" )
pairs( pcs[, 3:7] )
```

**Hint**: In Rstudio, you can click the 'zoom' button above above the plot to get a bigger version.

Our dataset contains samples with different ethnicities that in our dataset are recorded in the first column of the file (plink calls this a 'family ID' but here we've used it to record sampled ethnicities).  Make a list of them using table():

table( pcs[,1] )

Different ethnic groups might be genetically distinct, so they might separate on the PCs.  Let's replot colouring points by their ethnicity.  We'll do this in three steps

Turn ethnicities into integers that we can pass to the col argument of plot().
Plot using the colours
Add a legend to the plot so we can see what is what.

(The lower.panel and xpd arguments below are used to tweak the appearance of the plot and how the legend is plotted.  See ?pairs and ?par if you want more information, but don't worry about these for now.)

In RStudio:

```
pcs[,1] = as.factor(pcs[,1])
pcs$colour = as.integer(pcs[,1])
pairs(pcs[,3:7], col=pcs$colour, pch=20, lower.panel=NULL)
legend(0.1, 0.5, col=1:max(pcs$colour), legend = levels(pcs[,1]), pch=20, xpd=NA )
```

**Question**. Do the principal components reflect ethnicities in this dataset?  Which ethnicities are separated by which PCs? 

**Question**. Do any samples look especially genetically distinct from the other samples? What might be the reason for this?  Is the same set of samples outlying on all the PCs?  Use your R skills to identify these samples.  

# SNP weights/loadings

We would usually like to be able to interpret PCs as representing 'genome-wide' structure and ancestry of our samples. It's therefore wise to check the influence that each variant has on the principal components. To do this, let's examine the 
chr19-clean.eigenvec.var file that plink produced.  We'll load it into R and plot loadings across the chromosome.

In RStudio:
loadings = read.table("chr19-clean.eigenvec.var")
View(loadings)

Columns 3-22 of this file represent the loadings on PCs 1-20, respectively.  Let's plot loadings for the first 5 PCs.  (Here we use mfrow to make a plot with multiple rows.  The mar command adjusts the plot margins.  Again, see ?par for details, but don't worry about these for now.)

In RStudio:

```
par(mfrow=c(5,1), mar = c( 1, 4, 1, 2 ))
for( i in 1:5 ) {
plot( 1:nrow(loadings), abs(loadings[,i+2]), ylab=paste("PC ", i, sep="" ), ylim = c( 1, 10 ) )
}
```


**Question**. Do any SNPs stand out as having high loadings for particular PCs?  Which PCs?  Look back at your plot of principal components.  Do these PCs pick out particular clusters of samples?

**Question**. Go back and compute PCs without excluding related samples.  Then re-load and plot the PCs and loadings.  What drives the top PCs now?  What do the loadings look like?


# Plotting samples against a global dataset

A common strategy for identifying samples with unusual ancestry is to plot samples against a global reference dataset, such as that produced by the 1000 Genomes project (1000G).  A common way to do this is to compute a PCA of the external samples and project the GWAS dataset onto it (this can be done in plink using the `--pca` and `--pca-cluster-names` options - see the plink website for details).  For this practical, however, we'll take a simpler approach and compute a PCA of all the data together.  

To do this, we've created a merged dataset including our GWAS data with the African (AFR), European (EUR) and East Asian (EAS) samples from the 1000G reference panel.  To compute PCs,we sue a similar command as before, changing the dataset and exclusions files:

```
plink --vcf merged.with.1000G.vcf.gz --extract chr19-clean.prune.in --pca var-wts --out merged.with.1000G
```

As before this creates files names merged.with.1000G.eigenvec, etc.

Let's load the merged PCs and plot them with coloured ethnicities as before.  It's worth being careful with colours here, so we'll plot points with a particular colour scheme and use different shapes to distinguish the GWAS and reference panel samples.

In RStudio:

```
pcs = read.table( "merged.with.1000G.eigenvec" )
colnames(pcs)[3:7] = c("PC1","PC2","PC3","PC4","PC5")
pcs[,1] = factor(pcs[,1], levels = c("CAN","FAN", "JAN","RAN","AFR","EUR","EAS"))
palette = c("red2","green2","blue2","yellow3", "grey40","grey70","purple")
pcs$colour = palette[as.integer(pcs[,1])]
pcs$shape = 4
pcs$shape[pcs[,1] %in% c('AFR','EUR','EAS')] = 20
pairs(pcs[,3:7], col=pcs$colour, pch=pcs$shape, lower.panel=NULL)
legend(0.1, 0.5, col=palette, legend = levels(pcs[,1]), pch = c(rep(4,4),rep(20,3)), xpd=NA )
```7

**Question**. Which reference panel group do most of the GWAS dataset samples cluster with?  Are there any that don't?

**Question**. Use R to identify the samples that seem to cluster in the wrong place.  Do you recognise these samples?  Which reference panel group do they cluster near?  What do you conclude about these samples?

**Question**. Some of the 1000G samples labelled 'AFR'  also cluster nearer to the  Europeans than others.  Why is this?

**Hint**: Population codes can be found in the file resources/1000GP_Phase3.sample.  Look up the codes for a couple of samples, and look at the population definitions on the 1000 Genomes website at http://www.1000genomes.org/category/frequently-asked-questions/population.