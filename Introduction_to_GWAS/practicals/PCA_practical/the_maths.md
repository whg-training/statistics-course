[Up to the table of contents](Introduction.md) - [Back to the page on computing PCs](computing_PCs.md)

### An aside on the maths

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

#### Back to the computation...

When you're ready, [go back to the page on computing PCs](computing_PCs.md).

