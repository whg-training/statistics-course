[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_2.md)
[On to the next page](Getting_sequence_lengths.md)

Let's use our data to answer some questions on the extremes of the distribution of genes.

## Extremely big, small and complex genes

What is the largest gene, the gene with the most exons, or the gene with the most transcripts in each species?

You can work this out, but here are my solutions for reference:
[sql solution](solutions/part3/solutions.sql)
[python solution](solutions/part3/solutions.py)

Here are some things I noticed:

* In all the data I looked at, there's only one gene with > 100 transcripts. It is the human gene [*MAPK10*, "Mitogen-activated
  protein kinase 10"](https://www.uniprot.org/uniprot/P53779) and it has 151 transcripts! It even takes a while to load in
  [Ensembl](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000109339;r=4:85990007-86594625) and the
  [UCSC Genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A86016491%2D86594110&hgsid=275896231_HieWdPQTMOsgYQUFAnwTALgAECs0). If you explore these genome browsers you will find a wealth of
  information on this gene including its expression patterns and function. Is there much published literature about it?
  What does it look like in other species?

* Quite a few genes have lots of exons. One that really stands out in that [few species I've looked
  at](Counting_genes_2.md) is [*Titin*](https://en.wikipedia.org/wiki/Titin). It has 184 exons in Chimpanzees (*Pan
  troglodytes*) and ~113 in humans. Another is [*Nebulin*](https://en.wikipedia.org/wiki/Nebulin). These genes seem
  especially huge in Chimpanzees. Are they big in all species? In all great apes? Is there literature that might
  explain the results for these genes, or other genes with lots of exons?
  
* Those genes with lots of exons aren't the largest thoguh (by length). In humans and chimpanzees, that accolade goes
  to [`RBFOX1`]() with a whopping length of around 2 and a half megabases, or 0.08% of the genome. What does this look
  like across species? Where is it expressed? (Interestingly, this gene has been previously assocaited [with aggressive
  behaviour](https://www.nature.com/articles/s41380-018-0068-7)).

* There are also some pretty tiny genes around. In humans the smallest is
  [*ENSG00000288608*](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000288608;r=7:151061928-151061978
   ;t=ENST00000674552), which is only 50 base pairs long. (However, it seems to be contained in an exon of the much
  longer [*SLC4A2* Anion exhange protein
  2](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType
  =default&virtMode=0&nonVirtPosition=&position=chr7%3A151061877%2D151062029&hgsid=275886241_mxgWOGDr4elcf2SdW0zGz6ukwVJ
  W).  I'm not sure what this means!)

What other extreme genes are there?

### Single-exon genes

What proportion of protein-coding genes in each species have a single exon? (My solution is also in
[`solutions/part3/solutions.sql`](solutions/part3/solutions.sql) and
[`solutions/part3/solutions.py`](solutions/part3/solutions.py)). Almost 44% of [malaria
parasite](https://en.wikipedia.org/wiki/Plasmodium_falciparum) genes have a single exon, but the number is much lower
in other species I've looked at.

### Visualising

Visualisation (plotting) is beyond the scope of this tutorial - but you probably know how to do it already. It would be
interesting to visualise some of these statistics, for example

- histogram #transcripts, gene length, # of exons in each species
- plot the number of exons (or number of transcripts) versus gene length by species

and so on. The sky is the limit in terms of what you could investigate here - every gene has its own story. We'll talk
through your findings in the discussion session.

**Note.** Plotting is one thing you can't do in SQL - load the data into python or R for this. For reference, you can load
data into R using this type of code:

```
library( RSQLite )
db = dbConnect( dbDriver( "SQLite" ), "genes.sqlite" )
data = dbGetQuery( db, "SELECT * FROM gene_summary_view" )
```

### Getting sequence lengths

In the next section we will see how to [get sequence lengths](Getting_sequence_lengths.md).

