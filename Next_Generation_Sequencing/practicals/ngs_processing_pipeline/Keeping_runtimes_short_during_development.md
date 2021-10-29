[Click here to go back to the pipeline page](pipeline.md).

### Keeping runtimes short during pipeline development

When you're developing a pipeline, having a fast iteration time is important! You don't want to
wait two hours only to discover that it didn't work.

A good idea would therefore be to develop your practical using smaller, sub-sampled version of the datasets (whichever of the above raw data you use).  For example, you could run:

```
$ gunzip -c filename.fastq.gz | head -n 4000 | gzip -c > filename.subsampled.fastq
```

to take the first few reads from each file.

**Question.** The above command specifies a multiple of 4 lines.  Why?  How many reads does the above command extract?

If you set your pipeline up right (for example using a [config file](https://snakemake.readthedocs.io/en/stable/executing/cli.html) to list the data) then it will be easy to rerun it on the real data once you have it working.

