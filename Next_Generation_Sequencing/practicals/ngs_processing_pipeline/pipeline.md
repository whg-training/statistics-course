## Implementing a NGS data processing pipeline.

If you've followed the [introduction](README.md), you should now have a bunch of softare installed,
including the `snakemake` pipelining software, and you will have a set of 10 fastq files named in
the format `ERR[xxxxxx]_[1|2].fastq.gz`.  And you will also have downloaded the *P.falciparum* reference sequence, `Pf3D7_v3.fa.gz` - we recommend putting this in its own folder, say `data/reference/`.

**Challenge.**  Write a snakemake pipeline that processes these reads.  Your pipeline should:

* Take in a set of fastq read files named by accessions as described above.  To keep things well-organised, it's a good idea to keep these in a subdirectory, so they will look something like this:

```
   data/reads/ERR377582_1.fastq.gz
   data/reads/ERR377582_2.fastq.gz
   data/reads/ERR377591_1.fastq.gz
   data/reads/ERR377591_2.fastq.gz
   data/reads/ERR377629_1.fastq.gz
   data/reads/ERR377629_2.fastq.gz
   data/reads/ERR417621_1.fastq.gz
   data/reads/ERR417621_2.fastq.gz
   data/reads/ERR417627_1.fastq.gz
   data/reads/ERR417627_2.fastq.gz
```

* Output a set of BAM files containing these reads aligned to the a reference sequence.  The reads should be coordinate sorted, duplicates sohuld have been marked or removed, and the reads should be indexed.  My advice is to put results in a seperate subdirectory, so they will look something like this:

```
    results/aligned/QG0033-C.bam
    results/aligned/QG0033-C.bam.bai
    results/aligned/QG0041-C.bam
    results/aligned/QG0041-C.bam.bai
    results/aligned/QG0049-C.bam
    results/aligned/QG0049-C.bam.bai
    results/aligned/QG0056-C.bam
    results/aligned/QG0056-C.bam.bai
    results/aligned/QG0088-C.bam
    results/aligned/QG0088-C.bam.bai
```

* It should also start by QCing the fastq files. We suggest using `fastqc` and `multiqc` for this.
  So this will output some files that look like this:

```
    results/qc/QG0033-C.fastqc.html
    results/qc/QG0041-C.fastqc.html
    results/qc/QG0049-C.fastqc.html
    results/qc/QG0056-C.fastqc.html
    results/qc/QG0088-C.fastqc.html
    results/qc/multiqc_report.html
```
You should of course look at the output to look for anything odd!

* Optionally output a [`bedgraph`](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) file for each sample, reporting the coverage at each site in the genome.  (`bedtools genomecov -bg` is a good way to create this).  These will look like this:

```
    results/coverage/QG0033-C.coverage.bedgraph
    results/coverage/QG0041-C.coverage.bedgraph
    results/coverage/QG0049-C.coverage.bedgraph
    results/coverage/QG0056-C.coverage.bedgraph
    results/coverage/QG0088-C.coverage.bedgraph
```

* Optionally also output a variant calls file.  These will be an indexed bgzipped [vcf file](https://samtools.github.io/hts-specs/VCFv4.2.pdf), which should look like this:

```
    results/variant_calls/variant_calls.vcf.gz
    results/variant_calls/variant_calls.vcf.gz.tbi
```

To give you a better sense of this, here is a diagram of the overall pipeline you should implement:

![Diagram of pipeline](pipeline.svg).

Easy, right?

### Tips and tricks

Here is some guidance to help you write your pipeline.

#### How should I put sample information in?

I always put the sample information in through the config file - and specifically *not* encoded in
the snakefile itself. For example like this:

```
{
	"reference": "data/reference/Pf3D7_v3.fa.gz",
	"fastq_filename_template": "data/reads/subsampled/{id}_{read}.fastq.gz",
	"samples": [
		{ "name": "QG0033-C", "ID": "ERR377582" },
		{ "name": "QG0041-C", "ID": "ERR377591" },
		{ "name": "QG0049-C", "ID": "ERR417627" },
		{ "name": "QG0056-C", "ID": "ERR417621" },
		{ "name": "QG0088-C", "ID": "ERR377629" }
	]
}
```

And I then run snakemake with the `--configfile` argument like this:
```
snakemake -s pipelines/master.snakefile --configfile config.json
```

The point of this is that it makes it easy to run the pipeline on different sets of data


* [How should I put the sample information in?](How_to_get_sample_information_in.md)
* [How should I organise my pipeline files?](How_should_I_organise_my_pipeline.md)
* [My snakefiles are getting too big!](How_should_I_organise_my_pipeline.md)
* [Keeping runtimes short during development](Keeping_runtimes_short_during_development.md)
* [What's with all these intermediate files?](How_to_handle_intermediate_outputs.md)

### 

