[Click here to go back to the pipeline page](pipeline.md).

### My snakefiles are getting too big!

The [proposed pipeline](pipeline.svg) has quite a few steps, so your snakefiles probably will get pretty big.  A useful solution to this is as follows:

1. Split your snakefile into multiple snakefiles, one for each component.  (For example, you could have one for indexing the reference fasta file, one for aligning the data, one for variant calling, one for coverage computation, and so on.)

2. Use the snakemake ['include' statement](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html) to include them into a master snakefile.

For example, if you've structured your folder as described in [How should I organise my pipeline?](How_should_I_organise_my_pipeline.md), you might end up with this:

```
pipelines/
  master.snakemake
  reference.snakemake
  qc.snakemake
  alignment.snakemake
  variant_calling.snakemake
  coverage.snakemake
```

And your `master.snakefile` would look like this:

```
# contents of master.snakefile
include: "reference.snakemake"
include: "qc.snakemake"
include: "alignment.snakemake"
include: "variant_calling.snakemake"
include: "coverage.snakemake"

rule all:
  input:
    qc_report = "results/qc/fastq_multiqc_report.html",
    variants = "results/variants/variants.vcf.gz"
    # etc.
```

Then `qc.snakefile` would contain the rules for running `fastqc` and `multiqc`, `alignment.snakefile` would contain those for running the alignment steps, etc.

I find this quite useful as it makes it easy to drop bits of the pipeline in and out easily, and to keep each file relatively small and self-contained.

[Click here to go back to the pipeline page](pipeline.md).
