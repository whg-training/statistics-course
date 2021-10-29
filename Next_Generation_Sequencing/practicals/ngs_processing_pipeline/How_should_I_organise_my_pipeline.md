[Click here to go back to the pipeline page](pipeline.md).

### How should I organise my pipeline files?

It's worth organising your snakefiles right from the start. I typically use the following
conventions:

* All project data goes in the `data/` folder 
* All pipeline results go in the `results/` folder  (with minor exceptions e.g. some index files tend to get created next to their original files.)
* All pipeline files (snakefiles etc.) go in a seperate folder (I tend to call it `pipelines`).
* All the paths you use in the pipeline are expressed relative to the top-level folder.
* I then makes sure to always run the pipeline from the top-level folder.

This structure has several advantages:

* Because all the pipeline code is in a single place, it is easy to copy around / put up on github etc.
* All the results are also in one place, keeping things nice and clean.
* Because all the paths in the pipeline code are relative to the top-level directory, you don't have to worry about e.g. absolute pathnames etc.

Here's how it might look in practice for this tutorial:

```
  top-level folder/
    config.json

    pipelines/
      master.snakefile

    data/
      reads/
        ERR377582_1.fastq.gz
        ERR377582_2.fastq.gz
        (...)

      reference/
        Pf3D7_v3.fa.gz
        (...)

    results/
      (output files will go in here, in seperate folders)

    logs/
      (log files from the pipeline go here).
````

And I would run it like:
```
$ snakemake -s pipelines/master.snakefile --configfile config.json
```

**Note.** The snakemake documentation suggests a [similar, but slightly different layout](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

**Note.** If you find your snakefiles are getting too big, see [My snakefiles are too big!](My_snakefiles_are_too_big.md).
