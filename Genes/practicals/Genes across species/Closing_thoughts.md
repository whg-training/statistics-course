
[Go up to the table of contents](README.md)

## Closing thoughts

If you've made it through this tutorial - congratulations!  I hope it was useful to you.

If you haven't done any programming before then this will have been a bit of a challenge.  We've covered a number of topics - 

- basic programming in python
- using python data structures
- using the pandas data frame library
- writing command-line programs in python
- we've suggested some techniques for writing good, re-useable code
- using sqlite to store data
- writing simple algorithms like our [union-of-genomic-range algorithm](How_much_of_the_genome_is_in_genes.md)
- techniques for [dealing with memory issues](Memory_issues_and_how_to_solve_them.md) and for [scaling up analyses](Scaling_up.md) to larger datasets.

We haven't covered plotting (if you're a GMS student we'll have seperate tutorials about this) or working in other languages.

However, rather than just writing code, our analysis has been focussed on making sense of available gene annotation data.  For example we

- counted genes, transcripts and exons across species.
- got a sense of gene and genome complexity.
- measured genome coverage by protein-coding genes.
- recorded some data that starts to let us investigate specific genes in detail.

How far you take this  is up to you - there is a great deal of depth to this data and a lot remains to be discovered.

There is one **warning** I want to end on.  Not all genomes are assembled with the same accuracy, and not all genomes are annotated as well as humans.  Even for humans, some parts of the genome have not been assembled until recently - for example the [centromeres and telomeres](https://doi.org/10.1038/s41586-020-2547-7). We are analysing gene annotations not genes themselves.  To be really confident about comparisons you should also try to understand the larger context of your genome or genomic region of interest - and in particular how well annotated and assembled it is.  And there is lots of data out there to do it with...

