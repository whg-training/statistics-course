[Up to the table of contents](README.md) - [Back to the meta-analysis sectioin](meta-analysis.md))

### Using `FINEMAP` to fine-map associations

If you followed the [meta-analysis section](Meta-analysis.md) you should now have a dataframe with meta-analysis results.  It looks something like this:

    > meta_analysis
    # A tibble: 730 x 16
       rsid             chromosome position allele1 allele2 study1.beta study1.se   study1.P log10_study1.BF study2.beta study2.se study2.P meta.beta meta.se   meta.P log10_meta.BF
       <chr>                 <dbl>    <dbl> <chr>   <chr>         <dbl>     <dbl>      <dbl>           <dbl>       <dbl>     <dbl>    <dbl>     <dbl>   <dbl>    <dbl>         <dbl>
     1 19:49149607:G:A          19 49149607 G       A            0.202     0.228  0.376              -0.0925     -0.146      0.295  0.310      0.0176  0.0894 8.44e- 1       -0.0497
     2 19:49149711:G:A          19 49149711 G       A           -0.0183    0.202  0.928               0.0195      0.239      0.261  0.180      0.0245  0.0894 7.84e- 1       -0.0668
     3 19:49149937:T:G          19 49149937 T       G           -0.357     0.0785 0.00000540          4.90       -0.264      0.101  0.00463   -0.669   0.0894 7.45e-14       12.3   
     4 19:49150048:G:T          19 49150048 G       T           -0.332     0.0788 0.0000257           4.26       -0.271      0.102  0.00387   -0.636   0.0894 1.12e-12       11.2   
     5 19:49150130:G:T          19 49150130 G       T           -0.196     0.110  0.0756              0.947      -0.204      0.143  0.0767    -0.209   0.0894 1.96e- 2        1.52  
     6 19:49150233:G:A          19 49150233 G       A            0.133     0.134  0.321              -0.137      -0.0794     0.173  0.323      0.0381  0.0894 6.70e- 1       -0.0958
     7 19:49151312:C:T          19 49151312 C       T            0.238     0.145  0.101              -0.137       0.101      0.187  0.295      0.114   0.0894 2.03e- 1       -0.176 
     8 19:49151323:T:G          19 49151323 T       G            0.202     0.228  0.376              -0.0925     -0.146      0.295  0.310      0.0176  0.0894 8.44e- 1       -0.0497
     9 19:49151460:C:CA         19 49151460 C       CA          -0.283     0.153  0.0636              0.936      -0.218      0.197  0.134     -0.142   0.0894 1.12e- 1        0.830 
    10 19:49152118:C:T          19 49152118 C       T           -0.609     0.328  0.0632              0.634      -0.604      0.423  0.0767    -0.0723  0.0894 4.19e- 1        0.319 

We can look at the strongest signals in the data (either by P-value or Bayes factor) and compare across the studies and the meta-analysis.

However, **what if there are multiple causal variants in the region**?

If there are, then the single-SNP estimate at each SNP we see above would be a combination of several things:

* A contribution from the true, causal effect of the SNP.  (This will be zero if the SNP doesn't have a causal effect)
* A contribution from sampling noise (represented in the standard error)
* A contribution from other causal SNPs that happen to be in linkage disequilibrium (i.e. correlated) with this SNP.

Moreover it is possible that looking at the 'top' SNPs as above could be misleading in this case
(for example, if there were a SNP that is in strong LD with two causal SNPs, it might have an
apparent strong effect.)

In this part of the practical we will use the [`FINEMAP`](https:/www.finemap.me) tool to fit a model that allows for multiple possbile causal SNPs.


#### Preparing FINEMAP input files

To run finemap, we need to set up some input files in [the particular format that FINEMAP needs](https:/www.finemap.me).  `FINEMAP` takes:

* A set of per-SNP summary statistics (betas and standard errors) from the study. To maximise
  power, we will use your meta-analysed results here.

* A matrix of between-SNP linkage disquilibrium (correlation) values.

* It also requires us to specify the sample size - in our study the sample sizes were 12,500 (study
  1) and 7,500 (study 2), so the total sample size is 20,000.

We have precomputed the LD values for you in the file [`study_LD.ld`](study_LD.ld). (It is a giant
*730 &times; 730* matrix represented in a text file). But we need to make a file of meta-analysis
results.

If you look at [the documentation](https:/www.finemap.me) you'll see that FINEMAP needs a
space-separated file with these columns: rsid, chromosome, position, allele1, allele2, maf, beta,
se. The only thing we don't have that we need right now is the *minor allele frequency*. So let's
compute this first using the minor allele frequencies from the two studies and their sample sizes
(which in this data are study1: 12,500 and study2: 7,500):

```
meta_analysis$maf = ((study1$maf * 12500) + (study2$maf * 7500)) / (20000)
```

**Note.** The above calculation is **dangerous**: it would be wrong if the minor allele is not the
same in the two studies. For this data we happen to know that is not the case, but in real data you
should always aim to **compute the allele frequency of a specific allele** instead.

Now let's write a FINEMAP input file:

** In your R / RStudio session**:
```
finemap.data = meta_analysis[, c( "rsid", "chromosome", "position", "allele1", "allele2", "maf", "meta.beta", "meta.se" )]
colnames( finemap.data )[7:8] = c( "beta", "se" )
write_delim( finemap.data, "combined_study.z", delim = " ")
```

To get `FINEMAP` to run we also have to give it another input file telling it what all the input
and output files are. This one is a bit weird: it is a *semi-colon delimited* text file. Create a
new file called `finemap.master` and put this in:
```
z;ld;snp;config;cred;log;n_samples
combined_study.z;study_LD.ld;finemap.snp;finemap.config;finemap.cred;finemap.log;20000
```

**Note.** The format of the file is [described here](http://www.christianbenner.com/#input). What
input files are we using? What output files have we specified? What's the number `20000`?

Now you are ready to run `FINEMAP`

#### Running `FINEMAP``

I'm going to assuem you have downloaded and been able to run `FINEMAP` as described in the
[Introduction](Introduction.md). (If not please go and download it now.) Hopefully you know where
it is and can run it with a command of the following form. (**Note.** in this section we are back
to working on the command-line - these are not R commands!)

```
$ /path/to/finemap_v[version]/finemap_v[version]
```
On my older Mac OS system I had to use finemap v1.3, and running it looks like this:

    $ ~/Projects/Software/3rd_party/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX 

    |--------------------------------------|
    | Welcome to FINEMAP v1.3.1            |
    |                                      |
    | (c) 2015-2018 University of Helsinki |
    |                                      |
    | Help :                               |
    | - ./finemap --help                   |
    | - www.finemap.me                     |
    | - www.christianbenner.com            |
    |                                      |
    | Contact :                            |
    | - christian.benner@helsinki.fi       |
    | - matti.pirinen@helsinki.fi          |
    |--------------------------------------|
    
    Error : No flags detected!

For ease of use, let's copy the finemap executable into this directory:
```
cp /path/to/finemap_v[version]/finemap_v[version] ./finemap
```

Running finemap's is hopefully now easy:
```
$ ./finemap --sss --in-files finemap.master 
``` 

Be prepared for a bit of a wait as FINEMAP searches the space of possible causal configurations.

**Note.** What *is* FINEMAP doing? Well, you can read the [FINEMAP
paper](https://doi.org/10.1093/bioinformatics/btw018), but in short: finemap conducts a *shotgun
stochastic search* arond the space of possible causal configurations (groups of SNPs that might be
causal). It starts off with all single-SNP models, then it randomly tries a bunch of ways of adding
a second SNP to the first (picking the best model).  Then it keeps going randomly trying to add or
remove SNPs so that it explores the space of configurations that have high probability.

In general it would need the raw genotype data to do this. However, because it assumes an additive
model (of genotypes on phenotype) it turns out that an approximation based on the SNP effects and
ld is sufficient, which is how FINEMAP can run so fast.

#### Interpreting FINEMAP output

You run of FINEMAP should have produced several output files (the ones you named in your `finemap.master` file):

[TODO: improve this section].

* A 'configurations' file (`finemap.config`), listing the causal configurations that FINEMAP thinks are plausible.
* A 'credible intervals' file, listing variants wiht high probabilty of being the causal SNP for each signal in the top configuration.
* A 'SNP' file, listing the evidence that each SNP is one of the causal ones across configurations



[TODO] Fix the above / interpret here.










