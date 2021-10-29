## Required software

For the '**Introduction to NGS: from sequencing to variant calling**' workshop you will need the following software installed:

* For alignment: [`bwa`](https://github.com/lh3/bwa)
* For general data manipulation: [`samtools`](https://github.com/samtools/samtools)
* For BAM file processing: [`picard`](https://broadinstitute.github.io/picard/).

For the pipeline-building tutorial you will additionally need:

* For pipelining: [`snakemake`](https://snakemake.readthedocs.io/en/stable/).  (Or another workflow management tool of your choice.  [`WDL`](https://openwdl.org) and [`Nextflow`](https://www.nextflow.io) are possibilities, if you prefer them, but you will have to work them out for yourself.)
* For QC: [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`multiqc`](https://multiqc.info)
* For coverage calculations: [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html)
* For variant calling: [`octopus`](https://github.com/luntergroup/octopus)

With luck you should be able to install the above software using `conda` (please try now!) E.g. this may work better:

```
$ conda install samtools
```

```
$ conda install snakemake
```

and so on.

**Note.** I found that I generally needed to specify both `conda-forge` and `bioconda` channels to get these to work.  E.g.
```
$ conda install -c bioconda -c conda-forge snakemake fastqc multiqc samtools octopus bedtools
```

**Note.** There are also slightly different instructions on the [snakemake webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) that you might like to follow instead - this suggests to use the `mamba` command instead of `conda`.  (Mamba does seem slightly nicer - and it prints a picture of a snake!)

**Note.** Depending on your setup you might also need to use `sudo` to run these commands.  (This is related to system permissions - `sudo` gives a command higher priveleges.)  So you would write:
```
sudo conda install -c bioconda -c conda-forge snakemake
```

You will probably be asked for your system password at this point.

## Troubleshooting

We would like to identify and solve any installation issues.  If you are running this as part of a [WHG](https://www.well.ox.ac.uk) course please try this and report any problems back to us.

