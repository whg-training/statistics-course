[Up to the table of contents](Introduction.md) - [Back to the setup page](getting_setup.md) - [Forward to the page on LD pruning](ld_pruning.md)

### Overview of the practical

In this practical we will use `plink` to do several things to the data:

* to remove closely-related samples
* to compute principal components
* and to compute the SNP weights or loadings that tell us how principal components are weighted across the genome.

We'll also use `R` (we recommend [`RStudio`](https://www.rstudio.com)) to inspect and plot results.

Open a terminal window and first make sure you are in the right directory:

```
cd /path/to/PCA_practical
```

**Note.** This should be the folder where [you downloaded the data](getting_setup.md).

Also, in R / RStudio please set this directory as your current directory - either using the `setwd()` command like this:

*In R/RStudio:*
```
> setwd( '/path/to/PCA_practical' )
```

or by using the menu option `Session`->`Set working directory`->`Choose Directory`) in RStudio.

**Note.** The `>` indicates the `R` command prompt above - don't type that bit in!

### A note on quality control

Before carrying out a genetic analysis like PCA, it's important to have a good-quality dataset, and
this typically means carrying out careful quality control (QC) first. On this course we'll cover QC
in later lectures and practicals. For this practical we'll use an already-cleaned dataset contained
in the file `chr19-clean.vcf.gz`. You can look at the data in this file by typing

```
less -S chr19-clean.vcf.gz
```

**Note.** If you are using Mac OS X, you will need to use `zless` instead of `less` because the file is gzipped.

**Note.** Press `q` when you want to quit `less`.

The data consists of genotype calls at different sites (rows) for different samples (columns).  Feel free to look at the data by scrolling around using the arrow keys. When you've finished, press the 'q' key to quit back to the terminal prompt.

**Note.** This is a [Variant Call Format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file.
If you followed our earlier [next-generation sequencing
practical](../../../Next_Generation_Sequencing/practicals/ngs_processing_pipeline/) you will have
created a VCF file of variant calls from some sequence reads.  The VCF file for this practical is much simpler, though, because we have only included genotype calls.

### Preparing data for PCA

Before computing PCs we will need to do some pruning of the data.  We will:

* remove SNPs that are in high LD (to avoid confounding the analysis by local LD patterns.)
* remove samples that are too closely related (so that our PCs reflect the majority of our data.)

When you're ready, [go here to start pruning](ld_pruning.md).
