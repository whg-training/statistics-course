[Click here to go back to the pipeline page](pipeline.md).

### What's with all these intermediate output files?

All pipelines have a certain amount of *stuff* that gets generated along the way, but that you don't really want.  Read alignment pipelines are particularly bad for this - if you look at the [pipeline schema](pipeline.svg) you'll see that:

* the alignment step outputs a SAM file...
* which is converted to a BAM file...
* which you then have to sort by position...
* in which you then have to mark the duplicates...
* which are then indexed.

This only gets worse if you add in things like the adapter trimming step.

There are a couple of ways to handle this.  One way is to combine rules together and use UNIX pipes to send data from one command into the next.  For example:

```
bwa mem reference.fa read1.fq read2.fq | samtools fixmate -m - output.bam
```

That's ok, and in one sense it's a good idea since it reduces the number of temp files that have to
be made. But I find it can be annoying to debug when one of the commands fails (especially on a big input file), and it's not always possible.

The way that I sometimes like better is to make intermediate rules output files to a temporary location, and use the `temp()` function to tell snakemake we don't care about the intermediates.
This works well with snakemake's feature to *use the output of one rule as an input of another*.  Like this:

```
rule align_reads:
  input:
    fq1 = "data/read1.fq",
    fq2 = "data.read2.fq",
    ref = "data/reference/reference.fa"
  output:
    sam = temp( "results/alignment/tmp/aligned.sam" )
  shell: "bwa mem -o {output.sam} {input.ref} {input.fq1} {input.fq2}"

rule fix_matepair_information_and_convert_to_bam:
  input:
    sam = rules.align_reads.output.sam
  output:
    bam = temp( "results/alignment/tmp/aligned.bam" )
  shell: "samtools fixmate -m {input.sam} {output.bam}"
```

With this setup, all the temp files go into the `tmp/` subdir, but `snakemake` will realise they are temporary and delete them when they're no longer needed.

**Note.** A rule of thumb is that processing NGS data takes about 5 times as much space as the original data files - this is because of all these temporary files that have to be made.



