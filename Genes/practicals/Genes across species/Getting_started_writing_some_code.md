## Writing some code to process GFF files

Your task, if you choose to accept it, is to write some code that processes GFF files, make sense of them, and then
gathers statistics. You will then apply it to GFF files from multiple species to hopefully learn something about gene content of
the world's organisms.

We are planning to run this in the following way: during the course of this week, you work (on your own or in a group -
whatever works best) to develop some code to analyse a GFF file - and then to apply it to one or more GFF files to gain
some understanding of gene content across species. 

We're not writing this code for it's own sake but to answer our questions like the ones in our introduction:

- How many genes are there?
- How big are they?
- How much of the genome is in genes?
- How complex are genes - How many exons?  How many different transcripts?
- How much of genes is actually protein-coding sequence - and how much is untranslated?
- How much do these patterns differ across species?

Your code will have to parse the file and process it to extract meaning. And then you'll have figure out how to give
sensible quantitative answer to the above. We'll then discuss both the results and the code in the discussion session.

## Getting started

It might be that you already have a good idea how to go about this. If so, feel free to dive straight in. You're free
to use any language or system you like for this - standard options might be [python](https://www.python.org) or
[R](https://cran.r-project.org), but you could also use [julia](https://julialang.org), or even
[C++](https://en.wikipedia.org/wiki/C%2B%2B) or another compiled language. 

You are also free to use packages. However, for the rest of this tutorial we'll use [python](https://www.python.org) to
develop one way to solve this, without using any packages, and with an eye on our coding principles from [the
introduction](introduction.md).

