## Implementing a NGS data processing pipeline.

If you've followed the [introduction](README.md), you should now have a bunch of softare installed,
including the `snakemake` pipelining software, and you will have a set of 10 fastq files named in
the format `ERR[xxxxxx]_[1|2].fastq.gz`.

**Challenge.**  Write a snakemake pipeline that processes these reads.

Your pipeline should:

* take in a set of fastq read files named by accessions as above.  You should put these in a subdirectory, so they will look like this:

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

* output a set of BAM files containing these reads aligned to the a reference sequence.  The reads should be coordinate sorted, duplicates sohuld have been marked or removed, and the reads should be indexed.  My advice is to put results in a seperate subdirectory, so they will look something like this:

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

* also output a BED file for each sample, reporting the coverage at each site in the genome:

    results/coverage/QG0033-C.coverage.bed
    results/coverage/QG0041-C.coverage.bed
    results/coverage/QG0049-C.coverage.bed
    results/coverage/QG0056-C.coverage.bed
    results/coverage/QG0088-C.coverage.bed

* also (optionally) output a variant calls file.  These will come in a bgzipped vcf file, which should look like this:

    results/variant_calls/variants_calls.vcf.gz

* Along the way you will also need to QC the files - we suggest using `fastqc` and `multiqc`.  You might have something like this:

    results/qc/QG0033-C.fastqc.html
    results/qc/QG0041-C.fastqc.html
    results/qc/QG0049-C.fastqc.html
    results/qc/QG0056-C.fastqc.html
    results/qc/QG0088-C.fastqc.html
    results/qc/multiqc_report.html


Here are some of the things you will need in your pipeline:

* A QC step. We recommend running FastQC on the initial fastq files, and multiqc to get a report
  across the samples.

* A read alignment step.  For short-read paired-end data, `bwa mem` is the recommended algorithm

* A 'remove duplicates', or 'mark duplicates' step. This can either be done directly to the fastqc
  files (using `FastUniq`, which removes duplicate reads) or later at the alignment step (using
  `Picard MarkDup`, which doesn't remove them but flags them in the BAM file.)

* If your reads show high levels of adapter contamination - an adapter trimming step may be needed.
  For the current practical we suggest skipping this, but
  [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) is a good way to do this.

* A step to compute the coverage of your reads at each position in the genome - `bedtools genomecov` is a good option here.

Along the way you will need various other steps - sorting reads by alignment position, converting
between SAM and BAM formats, indexing and so on.

##

