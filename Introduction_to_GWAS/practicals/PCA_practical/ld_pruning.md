[Up to the table of contents](Introduction.md) - [Back to the previous page](overview.md) - [Forward to the next page](relatedness_pruning.md)

### LD pruning of SNPs

As described in the population genetics lecture this morning, genetic drift and other processes
lead to linkage disequilibrium (LD) between SNPs along a chromosome. To ensure the PCs we compute
represent genome-wide structure (not local LD) we'll first carry out LD pruning of our SNP set.

LD pruning removes correlated pairs of SNPs so that the remaining SNPs are roughly independent. (It
also helps to make subsequent computations quicker.) Run the following command to prune the dataset:

```
$ plink --vcf chr19-clean.vcf.gz --maf 0.01 --indep-pairwise 50 5 0.2 --out chr19-clean 
```

**Note**: commands like this should be typed out or pasted onto a single line in the terminal window (omitting the `$` of course!).  Then press <Enter> to run the command.  If all goes well you'll see a bunch of output on the screen, starting with some information about the version of plink that is running.

The above command tells plink to load the file `chr19-clean.vcf.gz` and to prune SNPs, leaving SNPs with minor allele frequency (MAF) at least 1%, and with no pairs remaining with pairwise r<sup>2</sup>>0.2.  (The other parameters, here 50 and 5, affect how the computation works in windows across the genome.  You can read about the behaviour here: [http://www.cog-genomics.org/plink2/ld](http://www.cog-genomics.org/plink2/ld)).

**Question**. Look at the screen output from the above plink command.  How many variants were in the original dataset?  How many were removed because their frequency was below 1%?  How many variants were removed due to LD pruning?  How many variants remain?

Type `ls` or use the file manager to view the directory.  The command above produced a number of files that all begin with the `chr19-clean` prefix.  For our purposes, the most important one is `chr19-clean.prune.in`, as this lists the SNPs that remain after pruning.  Feel free to look at all these files using less or a text editor.

#### Relatedness pruning
When you're ready, [go to the next page to identify and remove close relationships](relatedness_pruning.md).
