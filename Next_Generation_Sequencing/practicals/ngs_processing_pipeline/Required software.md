## Required software

For the 'Introduction to NGS: from sequencing to variant calling' workshop you will need the following software installed:

* For alignment: [`bwa`](https://github.com/lh3/bwa)
* For general data manipulation: [`samtools`](https://github.com/samtools/samtools)
* For BAM file processing: [`picard`](https://broadinstitute.github.io/picard/).

To put together a pipeline for the full tutorial you will additionally need:

* For pipelining: [`snakemake`](https://snakemake.readthedocs.io/en/stable/).  (Or another workflow management tool of your choice.  [`WDL`](https://openwdl.org) and [`Nextflow`](https://www.nextflow.io) are possibilities, if you prefer them, but you will have to work them out for yourself.)
* For QC: [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`multiqc`](https://multiqc.info)
* For coverage calculations: [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html)
* For variant calling: [`octopus`](https://github.com/luntergroup/octopus)

With luck you should be able to install the above software using `conda` (please try now!)  E.g.:

```
$ conda install snakemake
```

```
$ conda install fastqc multiqc
```

and so on.

Depending on your setup you might need to use `sudo` to run these commands.  (This is related to system permissions - `sudo` gives a command higher priveleges.)  So you would write:
```
sudo conda install snakemake
```

You will probably be asked for your system password at this point.

## Troubleshooting

We would like to identify and solve any installation issues.  If you are running this as part of a [WHG](https://www.well.ox.ac.uk) course please try this and report any problems back to us.
