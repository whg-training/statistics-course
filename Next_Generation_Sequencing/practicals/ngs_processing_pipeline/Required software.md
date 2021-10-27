## Required software

To run the tutorial you will need several pieces of software:

* For pipelining: [`snakemake`](https://snakemake.readthedocs.io/en/stable/).  (Or another workflow management tool of your choice.  [`WDL`](https://openwdl.org) and [`Nextflow`](https://www.nextflow.io) are possibilities, if you prefer them, but you will have to work them out for yourself.)
* For QC: [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`multiqc`](https://multiqc.info)
* For alignment: [`bwa`](https://github.com/lh3/bwa)
* For general file manipulation: [`samtools`](https://github.com/samtools/samtools)
* For variant calling: [`octopus`](https://github.com/luntergroup/octopus)

With luck you should be able to install these using `conda` (please try now!)  E.g.:

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
