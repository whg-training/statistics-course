[Up to table of contents](README.md)

[Back to the introduction](Introduction.md)

[Go the next page](Getting_started_writing_some_code.md)

## What gene annotation data looks like

The data we'll process is in a format called GFF. It has a specification that you can [read
here](https://m.ensembl.org/info/website/upload/gff3.html). 

Now is a good time to go and look at some of this data. You could download annotations from
[Ensembl](http://ftp.ensembl.org/pub/current_gff3/) for example. To make it easy to take a look,
I've pasted the top 1000 lines of a human gene annotation file in the file
`gencode.v38.annotation.head.gff`. This comes from the [GENCODE
project](https://www.gencodegenes.org)], which is a project that curates human and mouse genes).
And there are also the first 1000 lines of a *P.falciparum* gene annotation file in
`PlasmoDB-54_Pfalciparum3D7.head.gff` (from [PlasmoDB](https://plasmodb.org) in this folder. They
look a bit like this:

    ##gff-version 3
    #description: evidence-based annotation of the human genome (GRCh38), version 38 (Ensembl 104)
    #provider: GENCODE
    #contact: gencode-help@ebi.ac.uk
    #format: gff3
    #date: 2021-03-12
    ##sequence-region chr1 1 248956422
    chr1	HAVANA	gene	11869	14409	.	+	.	ID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2
    chr1	HAVANA	transcript	11869	14409	.	+	.	ID=ENST00000456328.2;Parent=ENSG00000223972.5;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-202;level=2;transcript_support_level=1;hgnc_id=HGNC:37102;tag=basic;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1
    chr1	HAVANA	exon	11869	12227	.	+	.	ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-202;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;hgnc_id=HGNC:37102;tag=basic;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1
    chr1	HAVANA	exon	12613	12721	.	+	.	ID=exon:ENST00000456328.2:2;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-202;exon_number=2;exon_id=ENSE00003582793.1;level=2;transcript_support_level=1;hgnc_id=HGNC:37102;tag=basic;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1

Like many bioinformatics file formats, this file doesn't bother to tell you what the columns are - or even what they
are called. You have to look at [the spec](https://m.ensembl.org/info/website/upload/gff3.html) for that. Some things
to figure out in there:

- How are genes represented?
- How are exons represented?
- How are transcripts represented?  (What is a transcript anyway?)
- How are exons, transcripts and genes linked in the file?

Can you work these out?

### What's in the file?

If you look at the spec you'll see the third column of the file above is called the `type`. The
spec says you'll see that this is part of a controlled 'ontology', which looks [a bit
scary](http://www.sequenceontology.org/so_wiki/index.php/Category:SO:SOFA). But what's actually in
the file? A quick way to find out is using the `cut` unix command, and the `sort -u` command to get
a unique list:

```sh
$ cut -f3 gencode.v38.annotation.head.gff | sort -u
$ cut -f3 PlasmoDB-54_Pfalciparum3D7.head.gff | sort -u
```

(Of course it would be much better to run this on the real files, not just the top 1000 lines I extracted above).

What you'll see in both files is that there is masses of information there! Including information on genes,
transcripts, exons, coding sequence (CDS), untranslated regions (UTR), translation start and stop codons.

**Question:** do the full files have any more types of record than these top-1000-lines files?  What are they?

If you want to see how this information shouldbe interpreted, try searching for a gene on
[Ensembl](http://www.ensembl.org/index.html) or the [UCSC Genome Browser](https://genome.ucsc.edu). E.g. for the gene
in the listing above it takes you to [this page on
ENSG00000223972](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000223972;r=1:11869-14409). This
already gives you a flavour of the complexity, because there are two transcripts on top of each other, with different
lengths and different numbers of exons.  

**Question:** Do the same for a protein-coding gene.  Does it look any different? Which website do you like better for
looking at genes, Ensembl or the [UCSC browser](https://genome.ucsc.edu)?  Spend some time exploring these sites.

## Write some code to process this

When you're ready, [write some code to process this data](Getting_started_writing_some_code.md).
