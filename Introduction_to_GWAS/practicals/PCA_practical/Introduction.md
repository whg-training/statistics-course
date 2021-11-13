# PCA practical

## Introduction

In this practical we will do principal component analysis (PCA), which is one of the fundamental tools of population genetics for identifying sample clustering and outliers. As explained in the talk by Gavin PCA is way of reducing a multidimensional space of data points to a set of orthogonal principal components that represent "directions" of greatest variance.  These may reflect population structure, but might also reflect cryptic relationships, poor QC, or other effects that lead to differences in genotype distribution.

### Getting set up

**Getting the software**. We will be using the software `PLINK` written by Christopher Chang:
[https://www.cog-genomics.org/](https://www.cog-genomics.org/)).  Before we start, please make sure you have downloaded this software and can run it in a terminal window on your system.  To check this, try running this in your terminal window:

```
$ plink
```

**Note.** The `$` in the above command indicates the command prompt - don't type that in!

You should see something like this:

    PLINK v1.90b6.17 64-bit (28 Apr 2020)          www.cog-genomics.org/plink/1.9/
    (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3

      plink <input flag(s)...> [command flag(s)...] [other flag(s)...]
      plink --help [flag name(s)...]

    Commands include --make-bed, --recode, --flip-scan, --merge-list, [...]

    "plink --help | more" describes all functions (warning: long).

**Getting the data**.  The data files for this practical can be found [in this folder](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Introduction_to_GWAS/practicals/PCA_practical/).  Please download these now.  For the first part of the practical we will use the `chr19-clean.vcf.gz` file, and later on we will use the other files as well.

For the practical we recommend making a new empty folder to put these in.  So when you're ready to go your folder should look something like this:

    PCA_practical/
        chr19-clean.vcf.gz
        merged.with.1000G.vcf.gz
        resources/
            1000GP_Phase3.sample

### The Practical

We will use `plink` to do several things:

* to remove closely-related samples
* to compute principal components
* and to compute the SNP weights or loadings that tell us how principal components are weighted across the genome.

We'll also use `R` (we recommend [`RStudio`](https://www.rstudio.com)) to inspect and plot results.

`plink`

To get started, open a terminal window and make sure you are in the right directory:

```
cd /path/to/PCA_practical
```

Also, in R / RStudio please set this directory as your current directory (either using the `setwd()` command or by using the menu option `Session`->`Set working directory`->`Choose Directory`).  Like this:

*In RStudio:*
```
> setwd( '/path/to/PCA_practical' )
```

**Note.** The `>` indicates the `R` command prompt here - don't type that bit in!

## A note on quality control

Before carrying out a genetic analysis like PCA, it's important to have a good-quality dataset, and
this typically means carrying out careful quality control (QC) first. On this course we'll cover QC
in later lectures and practicals. For this practical we'll use an already-cleaned dataset contained
in the file `chr19-clean.vcf.gz`. You can look at the data in this file by typing

```
less -S chr19-clean.vcf.gz
```

**Note.** If you are using Mac OS X, you will need to use `zless` instead of `less`.

The data consists of genotype calls at different sites (rows) for different samples (columns).  Feel free to look at the data by scrolling around using the arrow keys. When you've finished, press the 'q' key to quit back to the terminal prompt.

**Note.** This is a [Variant Call Format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file.
If you followed our earlier [next-generation sequencing
practical](../../../Next_Generation_Sequencing/practicals/ngs_processing_pipeline/) you will have
created a VCF file of variant calls from some sequence reads.  The VCF file for this practical is much simpler, though, because we have only included genotype calls.

## Preparing data for PCA

Before computing PCs we will need to do some pruning of the data.  We will:

* remove SNPs that are in high LD (to avoid confounding the analysis by local LD patterns.)
* remove samples that are too closely related (so that our PCs reflect the majority of our data.)

When you're ready, [go here to start pruning](ld_pruning.md).
