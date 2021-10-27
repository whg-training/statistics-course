## Next generation sequence data processing tutorial

In this tutorial / practical we are asking you to bring a set of [FASTQ
files](https://en.wikipedia.org/wiki/FASTQ_format) representing reads from an Illumina sequencing
platform through a basic sequence data processing pipeline. You will have to QC and align the
reads, identify or remove duplicate reads, compute coverage and - optionally generate an initial set of variant calls.

## Getting setup

### Required software.

You will need lots of software to implement this pipeline!  Before starting the tutorial please visit the [required software page](required software.md) and make sure you can install the needed pieces.

## Where to get the data

This tutorial will start with a set of fastq files representing *P.falciparum* (malaria) sequence reads from 5 samples.  The data comes from the [MalariaGEN Pf6 open resource](https://www.malariagen.net/resource/26) (described further in the [corresponding publication](https://wellcomeopenresearch.org/articles/6-42) and is publicly available via the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home).

To get the data you have two options:

1. Download the full FASTQ files from the ENA.  The relevant links are provided in [`samples.tsv`](samples.tsv).  The full data is aorund 12Gb.

2. To reduce the computational and space requirements for this tutorial, I created a lower depth version of the same data.  You can download these files from [this link](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/subsampled/).  the sub-sampled data is around 2Gb.

It's up to you which version of data you want to try to process.

### Some general advice

When you're developing a pipeline, having a fast iteration time is important! You don't want to
wait two hours only to discover that it didn't work.

A good idea would therefore be to develop your practical using smaller, sub-sampled version of the datasets (whichever of the above raw data you use).  For example, you could run:

```
$ gunzip -c filename.fastq.gz | head -n 4000 | gzip -c > filename.subsampled.fastq
```

to take the first few reads from each file.

**Question.** The above command specifies a How many reads does the above command take?

If you set your pipeline up right (for example using a [config file]() to list the data) then it will be easy to rerun it on the real data once you have it working.

