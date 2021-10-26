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

etc.  We would like to solve any issues so please try now!
