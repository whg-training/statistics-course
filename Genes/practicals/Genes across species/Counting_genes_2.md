[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_1.md)

## How protein-coding genes are there?

If you've followed so far you will have code `gff.py` that can parse a GFF file, and in the process
will pull out certain fields from the `attributes` column. This includes the `ID` attribute, the
`Parent` attribute that says how records are linked, and also the `biotype` column and maybe also
the `Name` (which contains gene names). In the [Ensembl
files](http://ftp.ensembl.org/pub/current_gff3/) the `biotype` is useful because it tells us what
kind of genes they are.

If you haven't already, please download the GFF files for some different species and run them
through your [sqlite conversion program](Converting_gff_to_sqlite.md) to get them into the
database. I'm going to assume you [followed the suggestion to add a column that distinguishes the
different species]() - in my code this is called `analysis` and I'll use that below.

(For reference my version of the code is at
[solutions/part2/gff_to_sqlite.py](solutions/part2/gff_to_sqlite.py) - feel free to run that
instead, if needed.)

This sqlite file is now in a good shape to start really exploring genes. First let's check out the
`biotypes`:

```
sqlite> .mode column
sqlite> .header on
sqlite> .width 50 25 25
sqlite> SELECT analysis, biotype, COUNT(*) FROM gff_data WHERE type=='gene' GROUP BY analysis, biotype ;
```
With the current data I have this gives:

    analysis                                            biotype                    COUNT(*)                 
    --------------------------------------------------  -------------------------  -------------------------
    Acanthochromis_polyacanthus.ASM210954v1.104         IG_J_gene                  2                        
    Acanthochromis_polyacanthus.ASM210954v1.104         IG_V_gene                  4                        
    Acanthochromis_polyacanthus.ASM210954v1.104         TR_J_gene                  5                        
    Acanthochromis_polyacanthus.ASM210954v1.104         protein_coding             24016                    
    Camelus_dromedarius.CamDro2.104.chr.gff3            IG_C_gene                  1                        
    Camelus_dromedarius.CamDro2.104.chr.gff3            IG_V_gene                  13                       
    Camelus_dromedarius.CamDro2.104.chr.gff3            TR_C_gene                  1                        
    Camelus_dromedarius.CamDro2.104.chr.gff3            TR_J_gene                  5                        
    Camelus_dromedarius.CamDro2.104.chr.gff3            TR_V_gene                  3                        
    Camelus_dromedarius.CamDro2.104.chr.gff3            protein_coding             18896                    
    Gallus_gallus.GRCg6a.104                            IG_V_gene                  98                       
    Gallus_gallus.GRCg6a.104                            protein_coding             16568                    
    Homo_sapiens.GRCh38.104                             IG_C_gene                  14                       
    Homo_sapiens.GRCh38.104                             IG_D_gene                  37                       
    Homo_sapiens.GRCh38.104                             IG_J_gene                  18                       
    Homo_sapiens.GRCh38.104                             IG_V_gene                  145                      
    Homo_sapiens.GRCh38.104                             TEC                        1056                     
    Homo_sapiens.GRCh38.104                             TR_C_gene                  6                        
    Homo_sapiens.GRCh38.104                             TR_D_gene                  4                        
    Homo_sapiens.GRCh38.104                             TR_J_gene                  79                       
    Homo_sapiens.GRCh38.104                             TR_V_gene                  106                      
    Homo_sapiens.GRCh38.104                             polymorphic_pseudogene     49                       
    Homo_sapiens.GRCh38.104                             protein_coding             19937                    
    Mus_musculus.GRCm39.104                             IG_C_gene                  13                       
    Mus_musculus.GRCm39.104                             IG_D_gene                  19                       
    Mus_musculus.GRCm39.104                             IG_J_gene                  14                       
    Mus_musculus.GRCm39.104                             IG_LV_gene                 4                        
    Mus_musculus.GRCm39.104                             IG_V_gene                  218                      
    Mus_musculus.GRCm39.104                             TEC                        3238                     
    Mus_musculus.GRCm39.104                             TR_C_gene                  8                        
    Mus_musculus.GRCm39.104                             TR_D_gene                  4                        
    Mus_musculus.GRCm39.104                             TR_J_gene                  70                       
    Mus_musculus.GRCm39.104                             TR_V_gene                  144                      
    Mus_musculus.GRCm39.104                             polymorphic_pseudogene     89                       
    Mus_musculus.GRCm39.104                             protein_coding             21834                    
    Pan_troglodytes.Pan_tro_3.0.104.chr                 IG_C_gene                  11                       
    Pan_troglodytes.Pan_tro_3.0.104.chr                 IG_V_gene                  91                       
    Pan_troglodytes.Pan_tro_3.0.104.chr                 TR_C_gene                  9                        
    Pan_troglodytes.Pan_tro_3.0.104.chr                 TR_V_gene                  66                       
    Pan_troglodytes.Pan_tro_3.0.104.chr                 protein_coding             21879                    
    PlasmoDB-54_Pfalciparum3D7                          protein_coding_gene        5318                     

The various biotypes used by Ensembl are [documented
here](https://m.ensembl.org/info/genome/genebuild/biotypes.html). As the above shows, the vast
majority of genes in the file are protein-coding genes. In each organism, however a few genes are
listed as `polymorphic_pseudogene`. You can [read more about pseudogenes
here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3491395/). They are genes that are non-coding in
the reference sequence (due to a disabling mutation), but still code for proteins in some
individuals.

The `IG_` and `TR_` biotypes are also interesting: they are the 'constant', 'joining', and
'variable' gene segments of immunuglobulin and T cell receptors. They do encode proteins, but via a
yet more [complex process that involves somatic recombination to assemble the mature genes in B and
T cells](https://en.wikipedia.org/wiki/V(D)J_recombination). These gene segments also lie in
regions that are [especially complex](https://doi.org/10.1371/journal.pcbi.1009254). 

(What are the `TEC` genes?)

However, the vast majority of the genes listed have type `protein_coding`, and we will focus on
these in the rest of this tutorial.

*Note.* Here is similar code to count genes across species in python:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM gff_data WHERE type IN ( 'gene' )", db )
data.groupby( [ 'analysis', 'biotype' ] ).size()
```

### What do some protein-coding genes look like?

**Challenge.** Pick a gene (say human ABO or FUT2) and investigate in detail using your file. (If you've followed so
far this will be in the `Name` column.) How many transcripts does it have? (Hint: find records with `Parent` equal to
the gene `ID`. (Are there any other records with the gene as parent, that aren't transcripts? What are they?). What is
the `transcript_support_level` for transcripts, and what does that mean? Pick a transcript and look at its sub-records
(again by matching the `Parent`` to the transcript `ID`). What types are they apart from exons? Can you see how this
data compares to the representation on [Ensembl](http://www.ensembl.org/) or on the [UCSC genome browser]()?

### What are all those species anyway?

If like me you're a bit unclear about what [all those
species](http://ftp.ensembl.org/pub/current_gff3/) are, now might be a good time to go and look at
the [OneZoom Tree of Life explorer](http://www.onezoom.org). This will tell you, for example, that
*Acanthochromis polyacanthus* is a member of the [Sticky Eggs
Group](http://www.onezoom.org/life/@Ovalentaria=5553750?img=best_any&anim=flight#x1307,y908,w1.5714)
 that also includes Cichlids, Silversides and Guppies, or that *mus musculus* is one of 37 species
collectively known as 'house mice'.

## How complicated are protein-coding genes?

Can we count how many transcripts each gene has, and how many exons?

**Challenge.** Write python code that reads in the appropriate data from `gff_data` table, and for
each gene counts i. the number of transcripts and ii. the average number of exons (averaged over
the transcripts for that gene).

**Hints.**

To do this is a bit more involved than previous steps. You need to find a way to join the
transcript records to the genes, and the exon records to the transcripts. Here is one way:
  
- Start by loading just the genes, transcripts, and exons into seperate data frames.

- Write a function `count_exons_per_transcript()` that takes the transcripts and exons, and
  returns a dataframe of transcripts with a column that reports the number of exons.
  
- Then similarly write a function `summarise_transcripts_per_gene()` that takes the genes and the
  above summarised transcripts, and returns a dataframe of transcripts with a column for the number
  of transcripts and a column for the average number of exons.

- There are a couple of ways to implement these functions. One way is to iterate through dataframes
  and build python data structures that capture the hierarchy of genes, transcripts, and exons.
  E.g. the first function might return something like this:
  
```
transcript_summary = [
   {
      "ID":"transcript:ENST00000641515",
      "exons": [
         "ENSE00003812156", "ENSE00003813641", "ENSE00003813949"
      ]
   },
   ...
}
```

and the second function something like this:

```
gene_summary = {
   "gene:ENSG00000186092": {
      "number_of_transcripts": 1,
      "average_number_of_exons": 3,
      "transcripts": [
         # array of objects recording transcripts and their exons
         {
            "ID":"transcript:ENST00000641515",
            "exons": [
               "ENSE00003812156", "ENSE00003813641", "ENSE00003813949"
            ]
         }
      ]
   },
   ...
}
```

You would then iterate through the second object to compute the summaries.

- A second way to implement this is to use 'data joins' and 'group by' operations. In pandas the
  functions to use are the [`merge()
  function`](https://pandas.pydata.org/docs/user_guide/merging.html)) and the
  [`groupby()`](https://pandas.pydata.org/docs/user_guide/groupby.html) function. A simple version
  of this is:

```
transcripts_and_exons = pandas.merge(
    transcripts,
    exons[ ["ID", "Parent" ]],
    how = "outer",
    left_on = "ID",
    right_on = "Parent"
)
transcript_summary = transcripts_and_exons.groupby( ['ID_x', 'Parent_x'] ).count()
```

This joins the transcripts to the exon IDs, and then groups the result by transcript to count the rows.

However the above has a few issues. First, gets the names wrong - "ID_x", "Parent_x" and so on, and
the third column does not even have a name! It is always good to have the right names. The
`.rename()` option and the 'Named aggregation' syntax described on the [`groupby()`
page](https://pandas.pydata.org/docs/user_guide/groupby.html) can be used to fix this. Something like this:

```
transcripts_and_exons.rename(
   columns = { 'ID_x': 'ID', 'Parent_x': 'Parent', 'ID_y': 'exon_ID', 'Parent_y', 'exon_Parent' },
   inplace = True
)

def count( x ):
   return len(x)

transcript_summary = (
   transcripts_and_exons
      .groupby( ['ID', 'Parent'] )
      .agg(
         number_of_exons = pandas.NamedAgg(
            column = "exon_ID",
            aggfunc = count
         )
      )
)
```

If you run this you will see it has "ID", "Parent" and "number_of_exons" columns. 

More seriously, the above code has a possible bug. (You did test it, right?) Specifically it gets
the answer wrong if a transcript has no exons. How did I discover this? [By writing a
test](solutions/part2/test_gff.py).  So that needs to be fixed too.  Good luck! 

**Note.** My version of the code can be found in the [`solutions/part2/gff.py`](solutions/part2/gff.py) file. There is
both a pandas version (`summarise_genes()`) and a python datastructure version (`summarise_genes_python_version()`).

**Note.** which of these approaches do you find easier to understand? The python version of my code is definitely
longer, but neither seems especially simple. However, the `count_exons_per_transcript()` and
`summarise_transcripts_per_gene()` functions are pretty similar and so they are good candidates for a refactor - so I
think this code can be improved...

### A sqlite approach

In the rest of this section I'll show how you could solve those problems in the sqlite database
itself using data joins and group by operations.

First, it's convenient to make some views of the data that just reflect the genes, transcripts and
exons:

```
CREATE VIEW gene_view AS SELECT * FROM gff_data WHERE type == 'gene' ;
CREATE VIEW transcript_view AS SELECT * FROM gff_data WHERE type == 'mRNA' ;
CREATE VIEW exon_view AS SELECT * FROM gff_data WHERE type == 'exon' ;
```

Second, to make lookups efficient, we'll also need to make sure the `ID` field is indexed:

```
CREATE INDEX gff_data_ID_index ON gff_data( ID ) ;
```
(This can take a minute or so as it builds the index).

Now we can join up the data to make the counts.  First we'll count exons in transcripts:
```
CREATE VIEW transcript_summary_view AS
SELECT T.*, COUNT(E.Parent) AS number_of_exons
FROM transcript_view T
LEFT JOIN exon_view E ON E.Parent == T.ID
GROUP BY T.ID;
```


If you look at the result:
```
SELECT * FROM transcript_summary_view LIMIT 10 ;
```
you'll see it takes a little while to run its computation, and it has one row per transcript, with
the number of exons in the last column.

**Note.** As with python in principle we have to be careful what we count above. The `COUNT(E.Parent)` makes sure that
we count only rows that correspond to an actual exon (rather than a transcript with no exons, which should get a count
of zero. (There shouldn't really be any transcripts like this - otherwise it isn't really a transcript - but it is
worth making sure. Are there any?)

Next we can summarise transcripts for each gene:
```
CREATE VIEW gene_summary_view AS
SELECT
    gene_view.*,
    COUNT(T.Parent) AS number_of_transcripts,
    CAST( SUM(number_of_exons) AS FLOAT ) / COUNT(T.ID) AS average_number_of_exons
FROM gene_view
LEFT JOIN transcript_summary_view T ON T.Parent == gene_view.ID
GROUP BY gene_view.ID;
```

This solves our problem:

```
SELECT ID, Name, biotype, number_of_transcripts, average_number_of_exons  FROM gene_summary_view LIMIT 10;
```

**Note.** Using `.mode line` or `.mode column` and `.header on` is useful in sqlite3. You can also
just run these queries and load them into python as described above.

Unfortunately these queries are taking... quite a while. In fact, we are asking sqlite to do quite
a bit of work here, largely because we've built this on dynamic views of the original `gff_data`
table - which is very big and includes plenty that is not relevant to our query. To speed up our
queries (in a way that mimics what we would do in python) we can 'bake' the relevant view into
their own tables first:

```
CREATE TABLE genes AS SELECT analysis, ID, Parent, Name, biotype, seqid, start, end, strand FROM gene_view ;
CREATE TABLE transcript_summary AS SELECT analysis, ID, Parent, Name, biotype, seqid, start, end, strand, number_of_exons FROM transcript_summary_view ;
CREATE INDEX transcript_summary_Parent_INDEX ON transcript_summary( Parent ) ;
```

and then rewrite `gene_summary_view` to use these tables:
```
DROP VIEW gene_summary_view ;
CREATE VIEW gene_summary_view AS
SELECT
    genes.*,
    COUNT(T.ID) AS number_of_transcripts,
    CAST( SUM(number_of_exons) AS FLOAT ) / COUNT(T.ID) AS average_number_of_exons
FROM genes
LEFT JOIN transcript_summary T ON T.Parent == genes.ID
GROUP BY genes.ID;
```
On my system, with the data above this brings query time down to a few seconds.

**Note.** The `CAST( .. AS FLOAT )` syntax is needed above to make sure sqlite uses [floating-point arithmetic]() for the division.

Let's use this data to [investigate genes more deeply](Counting_genes_3.md).

