[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_1.md)

## How many protein-coding genes are there?

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

### What are all those species anyway?

If like me you're a bit unclear about what [all those
species](http://ftp.ensembl.org/pub/current_gff3/) are, now might be a good time to go and look at
the [OneZoom Tree of Life explorer](http://www.onezoom.org). This will tell you, for example, that
*Acanthochromis polyacanthus* is a member of the [Sticky Eggs
Group](http://www.onezoom.org/life/@Ovalentaria=5553750?img=best_any&anim=flight#x1307,y908,w1.5714)
 that also includes Cichlids, Silversides and Guppies, or that *mus musculus* is one of 37 species
collectively known as 'house mice'.

## How complicated are genes?

Can we count how many transcripts each gene has, and how many exons?

**Challenge.** Write python code that reads in the appropriate data from `gff_data` table, and for
each gene counts i. the number of transcripts and ii. the average number of exons (averaged over
the transcripts for that gene).

**Hints.**

- To do this is a bit more involved than previous steps. We need to join the transcripts to the
  genes and the exons to the transcripts somehow. There are several ways to do this.

- In what follows I'll give a SQL approach to this - and then show how to mimic that approach in
  python. This uses data joins (implemented in the pandas [`merge()
  function`](https://pandas.pydata.org/docs/user_guide/merging.html)) to join exons to transcripts
  and transcripts to genes.
  
- But you don't have to do it this way. For example, another (better?) way would be to iterate through the data
  (for example using `.apply()`), and directly build a data structure mapping genes to transcripts
  and transcripts to exons. E.g. you might end up with a structure like this:

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

(You could also save quite a bit of it by storing row indices instead of gene and transcripts
identifiers in the above.)

### A sqlite approach

Just for completeness, here's how you could solve that in the sqlite database itself.

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
SELECT T.*, COUNT(*) AS number_of_exons
FROM transcript_view T
LEFT JOIN exon_view E ON E.Parent == T.ID
GROUP BY T.ID;
```

If you look at the result:
```
SELECT * FROM transcript_summary_view LIMIT 10 ;
```
you'll see it takes a little while to run its computation, and it has one row per transcript, with
the number of exons in the last column

Next we can summarise transcripts for each gene:
```
CREATE VIEW gene_summary_view AS
SELECT G.*, COUNT(*) AS number_of_transcripts, SUM(number_of_exons) / (COUNT(*) + 0.0) AS average_number_of_exons
FROM gene_view G
LEFT JOIN transcript_summary_view T ON T.Parent == G.ID
GROUP BY G.ID;
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
SELECT G.*, COUNT(*) AS number_of_transcripts, (0.0+SUM(number_of_exons))/COUNT(*) AS average_number_of_exons
FROM genes G
LEFT JOIN transcript_summary T ON T.Parent == G.ID
GROUP BY G.ID;
```
On my system, with the data above this brings query time down to a few seconds.

Equivalent python code [is here](solutions/part2/gene_summary.py).

Let's use this to answer some questions.

### Which genes have the most transcripts?

I'll use the SQL version to answer these questions.

**Challenge.** replicate these queries using python.

```

SELECT analysis, MAX(number_of_transcripts), MAX( average_number_of_exons ) FROM gene_summary_view
GROUP BY analysis ;

```

It turns out that (among the above data) humans have by far the largest number of annotated
transcripts for a single gene - 151, for the MAPK10 gene which is a [member of the MAP kinase
family](https://en.wikipedia.org/wiki/MAPK10). However, Chimpanzees turn out to have a gene with an
extremely large number of exons: 184 for the *TTN* gene which encodes
(Titin)[https://en.wikipedia.org/wiki/Titin]. (*TTN* also has the highest number of exons in
humans - but only a paltry 112.)  These observations may be explained by the description of Titin:

    Titin /ˈtaɪtɪn/, also known as connectin, is a protein that in humans is encoded by the TTN
    gene. Titin is a giant protein, greater than 1 µm in length,that functions as a molecular
    spring which is responsible for the passive elasticity of muscle. It comprises 244 individually
    folded protein domains connected by unstructured peptide sequences. These domains unfold when
    the protein is stretched and refold when the tension is removed.

    - Wikipedia

### Joining data in python

Here is the same process carried out in python.
First we load the data:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
genes = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'gene' )", db )
transcripts = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'mRNA' )", db )
exons = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand FROM gff_data WHERE type IN ( 'exon' )", db )
```

Now we create the per-transcript summary of the number of exons. For this we have to use the slightly complicated
pandas syntax for [`groupby` operations with named columns](https://pandas.pydata.org/docs/user_guide/groupby.html):

```
transcript_summary = pandas.merge(
   transcripts,
   exons[["ID", "Parent"]], 
   how = "outer",
   left_on = "ID",
   right_on = "Parent"
).groupby( ['ID_x', 'Parent_x'], as_index = False ).agg(
    number_of_exons = pandas.NamedAgg( column = "ID_y", aggfunc = numpy.size )
)
```

And then we create the per-gene summary:

```
gene_summary = pandas.merge(
    genes,
    transcript_summary,
    how = "outer",
    left_on = "ID",
    right_on = "Parent_x"
).groupby( [ "analysis", "ID" ] ).agg(
    number_of_transcripts = pandas.NamedAgg( column = "ID", aggfunc = numpy.size ),
    average_number_of_exons = pandas.NamedAgg( column = "number_of_exons", aggfunc = numpy.mean ),
)

gene_summary
```

This prints:

                                                                         number_of_transcripts  average_number_of_exons
    analysis                                    ID                                                                     
    Acanthochromis_polyacanthus.ASM210954v1.104 gene:ENSAPOG00000000002                      1                      6.0
                                                gene:ENSAPOG00000000003                      1                      6.0
                                                gene:ENSAPOG00000000004                      1                     10.0
                                                gene:ENSAPOG00000000005                      1                      9.0
                                                gene:ENSAPOG00000000006                      1                      3.0
    ...                                                                                    ...                      ...
    PlasmoDB-54_Pfalciparum3D7                  PF3D7_API04600                               1                      1.0
                                                PF3D7_API04700                               1                      1.0
                                                PF3D7_MIT01400                               1                      1.0
                                                PF3D7_MIT02100                               1                      1.0
                                                PF3D7_MIT02300                               1                      1.0
    
    [134092 rows x 2 columns]

You can also check this gives the same results as the above SQL:

```
gene_summary_by_sql.set_index(['analysis', 'ID'], inplace = True )
gene_summary.sort_index( inplace = True )
gene_summary_by_sql.sort_index( inplace = True )

comparison = gene_summary[['number_of_transcripts', 'average_number_of_exons']] == gene_summary_by_sql[['number_of_transcripts', 'average_number_of_exons']]
```

## Sequence lengths

