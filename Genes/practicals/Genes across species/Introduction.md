[Up to table of contents](README.md)

# Coding warm-up - analysing genes across species.

One of the most amazing scientific acheivements in the last 100 years has been the mapping of genes - DNA
sequences that encode proteins or other functional molecules - across many hundreds of species. For example, the
[Ensembl ftp site](http://ftp.ensembl.org/pub/current_gff3/) holds gene annotation files for over 300 species, from
[*Acanthochromis Polyacanthus*](https://en.wikipedia.org/wiki/Spiny_chromis) to [*Zosterops
Lateralis*](https://en.wikipedia.org/wiki/Silvereye), while other sites hold genes for other organisms such as the [*Plasmodium*
parasites](https://plasmodb.org/plasmo/app/downloads/Current_Release/) that cause malaria.

But what do genes actually look like?  Lots of questions spring to mind:

- How many genes are there?
- How big are they?
- How much of the genome is in genes?
- How complex are genes - How many exons?  How many different transcripts?
- How much of genes are actually protein-coding - and how much is untranslated?
- Do these patterns differ across species?  How?

You can probably think many others!

In this practical we will write some re-useable code that can process gene annotation files, as a way to start figuring
out these questions.

### Plan of attack

This problem is fairly typical of bioinformatics problems in the following way. There is some data
available (like those in the above FTP sites). It comes in file formats that somebody has invented.
And we have a set of slightly vague scientific questions in mind we're interested in. To answer
them, we have to understand the data files, write some code to process it, and come up with some
sort of quantitative analysis.

## Coding for success

In this tutorial we'll focus on two things at once. Our primary focus will be on answering the
scientific questions above. But we're also going to use this as a chance to write some useful and
re-useable code.

There are lots of ways to define 'good code' and lots of ways to write it. But here are some simple
things our code ought to aim for:

- it ought to work
- it ought to not take too long to do it
- it ought to be obvious what it does

Arguably the third point is the most important one here. Because if you can figure out what the code does, then you
can fix any problems with it and work to speed it up. But if you can't figure it
out, you'll end up throwing it away and starting again.

Here are some useful principles to keep in mind as we try to write useful code:

1. *Keep it simple.*
2. *Keep related things together.*
3. *Don't repeat yourself.*
4. *Write testable code.*
5. *Write tests.*
6. *Write short functions.*
7. *Give things good names.*

There are of course lots of other principles that we could apply, but these are ones I tend to focus on.

### A comment on comments

Many people would add to the above: 

7. *Write lots of comments to explain the code.*

However, I don't generally do this. This is because I think if 1-6 above are done correctly, then the
code will largely explain itself (at least to a trained eye). I
try to reserve comments for the inevitable leftover bits of code that are still surprising or hard to
understand, even after I've done my best on the above. (Not writing comments also has the great advantage of being
less work, especially when you start to change the code... have you remembered to change the comment?)

## What gene annotation data look like

When you've had enough coding philosophy, [go and see what gene annotation data looks like](What_gene_annotation_data_looks_like.md).
