[Up to the table of contents](Introduction.md) - [Back to the previous page](computing_PCs.md) - [Forward to the next page](plotting_against_a_global_dataset.md).

### Plotting samples against a global dataset

A common strategy for identifying samples with unusual ancestry is to plot samples against a global
reference dataset, such as that produced by the 1000 Genomes project (1000G). A common way to do
this is to compute a PCA of the external samples and project the GWAS dataset onto it. In plink,
this can be done in plink using the `--pca` and `--pca-cluster-names` options - see the [plink
website](https://www.cog-genomics.org/plink/) for details). For this practical, however, we'll take
a simpler approach and compute a PCA of all the data together.

To do this, we've created a merged dataset `merged.with.1000G.vcf.gz` including our GWAS data with the African (AFR), European (EUR) and East Asian (EAS) samples from the 1000G reference panel.  To compute PCs, we use a similar command as before, changing the dataset and exclusions files:

```
plink --vcf merged.with.1000G.vcf.gz --extract chr19-clean.prune.in --pca var-wts --out merged.with.1000G
```

As before this creates files names `merged.with.1000G.eigenvec`, etc.

Let's load the merged PCs and plot them with coloured ethnicities as before.  It's worth being careful with colours here, so we'll plot points with a particular colour scheme and use different shapes to distinguish the GWAS and reference panel samples:

In R/RStudio:

```
pcs = read.table( "merged.with.1000G.eigenvec" )
colnames(pcs)[1:7] = c( "group", "ID", "PC1", "PC2", "PC3", "PC4", "PC5" )

palette = c(
  CAN = "red2", FAN = "green2", JAN = "blue2", RAN = "yellow3",
  AFR = "grey40", EUR = "grey70", EAS = "purple"
)
shapes = c(
  CAN = 4, FAN = 4, JAN = 4, RAN = 4,
  AFR = 20, EUR = 20, EAS = 20
)

pcs$colour = palette[pcs$group]
pcs$shape = shapes[pcs$group]

pairs( pcs[,3:7], col = pcs$colour, pch = pcs$shape, lower.panel = NULL )
legend(
  0.1, 0.5,
  col = palette,
  legend = names(palette),
  pch = shapes,
  xpd=NA
)
```

**Question**. Which reference panel group do most of the GWAS dataset samples cluster with?  Are there any that don't?

**Question**. Use R to identify the samples that seem to cluster in the wrong place.  Do you recognise these samples?  Which reference panel group do they cluster near?  What do you conclude about these samples?

**Question**. Some of the 1000G samples labelled 'AFR' also cluster nearer to the Europeans than others.  Why is this?

**Hint**: Population codes can be found in the file `resources/1000GP_Phase3.sample`.  Look up the codes for a couple of samples, and look at the population definitions on the 1000 Genomes website at `http://www.1000genomes.org/category/frequently-asked-questions/population`.

#### Summary

Congratulations!  You have completed a basic principal components analysis.  We will use these PCs later when we run the [genome-wide association study practical]("../GWAS_analysis_practical").


