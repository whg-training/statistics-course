[Click here to go back to the pipeline page](pipeline.md).

### How should I get sample information into the pipeline?

You have sample information in the files [`samples.tsv`](samples.tsv). You could hard-code that directly into the snakefile, like this:

```
# master.snakefile

samples= [
	{ "name": "QG0033-C", "accession": "ERR377582" },
	{ "name": "QG0041-C", "accession": "ERR377591" },
	{ "name": "QG0049-C", "accession": "ERR417627" },
	{ "name": "QG0056-C", "accession": "ERR417621" },
	{ "name": "QG0088-C", "accession": "ERR377629" }
]

# rules go here
```

```
import pandas
samples = pandas.read_table( "samples.tsv" )
```

However these options are annoying because they hard-code the data into your snakefile.  If you want to run the pipeline on some other data, or a test datasets, or a larger or smaller dataset, you're stuck with editing the file.

Instead, what I do is use *config file* to input the data.  For example, for this project my config file looks like this:
```
{
	"reference": "data/reference/Pf3D7_v3.fa.gz",
	"fastq_filename_template": "data/reads/{accession}_{read}.fastq.gz",
	"samples": [
		{ "name": "QG0033-C", "accession": "ERR377582" },
		{ "name": "QG0041-C", "accession": "ERR377591" },
		{ "name": "QG0049-C", "accession": "ERR417627" },
		{ "name": "QG0056-C", "accession": "ERR417621" },
		{ "name": "QG0088-C", "accession": "ERR377629" }
	]
}
```

I then run snakemake passing in this config file, like this:
```
$ snakemake -s pipelines/master.snakefile --configfile config.json
```

This makes it easy to do a test run: just use a second, test config file (`config_test.json`, say).

**Note.** This works well for this handful of samples.  If there were lots of samples I would probably instead put the name of the `samples.tsv` file into the config file.

