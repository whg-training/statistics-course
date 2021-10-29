## Next generation sequence data processing tutorial

In this tutorial / practical we are asking you to bring a set of [FASTQ
files](https://en.wikipedia.org/wiki/FASTQ_format) representing reads from an Illumina sequencing
platform through a basic sequence data processing pipeline. You will have to QC and align the
reads, identify or remove duplicate reads, compute coverage and - optionally generate an initial set of variant calls.

## Getting setup

### Required software.

You will need lots of software to implement this pipeline!  Before starting the tutorial please visit the [required software page](required software.md) and make sure you can install the needed pieces.

### Getting the data

This tutorial will start with a set of fastq files representing *P.falciparum* (malaria) sequence reads from 5 samples.  The data comes from the [MalariaGEN Pf6 open resource](https://www.malariagen.net/resource/26) (described further in the [corresponding publication](https://wellcomeopenresearch.org/articles/6-42) and is publicly available via the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home).

To get the data you have two options:

1. Download the full FASTQ files from the ENA.  The relevant links are provided in [`samples.tsv`](samples.tsv).  The full data is around 12Gb.

2. To reduce the computational and space requirements for this tutorial, I created a lower depth version of the same data.  You can download these files from [this link](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reads/subsampled/).  The sub-sampled data is around 2Gb.

It's up to you which version of data you want to try to process.

You will also need the reference sequence to align to - we will use `Pf3D7_v3` which you can find [here](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reference/).

## Getting started

When you're ready, go the the [pipeline overview page](pipeline.md).