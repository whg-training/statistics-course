[Up to table of contents](Introduction.md)
[Back to the previous page](Converting_gff_to_sqlite.md)

## Gathering some first statistics

Let's start to attack our scientific questions. From simplicty, from here on in we will work with the GFF3 files
from the [Ensembl ftp site](http://ftp.ensembl.org/pub/current_gff3/) - also those from
[PlasmoDB](https://plasmodb.org/plasmo/app/downloads/Current_Release/). These use the terminology `mRNA` for a gene
transcript, and they also have the genome sequence lengths written in the metadata, making life easy for us.

We will focus on protein-coding genes, and their transcripts, exons and coding sequence. They have type `gene`, `mRNA`, `exon`
and `CDS` in the files respectively. They come in a basic hierarchy:

- Each gene can have multiple transcripts (i.e. multiple expressed forms at mRNA level - e.g. they might differ in how
  exons are spliced together, or have different transcription start or end sites).
  
- Each transcript is made up of one or more exons. The introns in between are spliced out. (This happens during
  transcription process.)

- Of the exon sequence that makes it into a transcript, only a portion is actually translated into the mature protein
  during processing by the ribosomes.  This is the *coding sequence* and the rest is *untranslated sequence*.

The files we're looking are (roughly speaking) humanities' best guess at what this picture looks like in each organism.


We're actually in a good shape for many of [our scientific questions](Introduction.md).  Let's take them one by one:

### How many genes are there?

Let's start with counting genes - easy!  One way is to use sqlite.  For example using the sqlite file you created in the last step:

```
sqlite> SELECT COUNT(*) FROM genes WHERE type==`gene` ;
```

*Note.* I ran the above command interactively in sqlite3 - you get there by typing `sqlite3 genes.sqlite` in the shell.
Type `.quit` to quit.)

If you followed the suggestions, or used my version fo the code, your file will also have an `analysis` column that
lets you differentiate species.  So you can make these counts for multiple species:

```
sqlite> .mode column
sqlite> .header on
sqlite> SELECT analysis, COUNT(*) FROM genes WHERE type=='gene' GROUP BY analysis ;
```

Which gives:

    analysis                                     COUNT(*)  
    -------------------------------------------  ----------
    Acanthochromis_polyacanthus.ASM210954v1.104  48054     
    Homo_sapiens.GRCh38.104                      21451     
    Mus_musculus.GRCm39.104.chr.gff3             25655     


You can also do this pretty easily in python:
```
>>> import gff
>>> data = gff.parse_gff3_to_dataframe("Homo_sapiens.GRCh38.104.chr.gff3" ) # assuming you have downloaded this file
>>> data[ data['type'] == 'gene' ].shape
(21451, 11)
```

**Warning.** When I do this the process is using up 1Gb of memory (as reported by `top`).  That's second only to Microsoft Word!

**Question:** how could you best replicate the above across-species query in python?

The above counts all "genes" - not all of these code for proteins though.  

**Question:** How do we count protein-coding genes?

**Challenge:** Can you update your `gff.py` code to make it easy to get this information out?

**Hints:** 

- There's a lot of stuff in the `attributes` column - what looks useful?
- You've already extracted `ID` and `Parent` into dedicated columns - it ought to be easy to get other attributes as well?
- (I'd suggest getting the gene `Name` as well where you're about it.)

**More hints:** My solution is in [`solutions/part2`](solutions/part2/) if you want to take a look.

**Challenge:** make a table of counts of genes by type across several species.

My attempt at this is described [on the next page](Counting_genes_2.md).
