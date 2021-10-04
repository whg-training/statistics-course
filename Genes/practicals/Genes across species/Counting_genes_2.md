[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_1.md)

## How many protein-coding genes are there?

If you've followed so far you will have code `gff.py` that can parse a GFF file, and in the process
will pull out certain fields from the `attributes` column. This includes the `ID` attribute, the
`Parent` attribute that says how records are linked, and maybe also the `Name` column (gene names)
and the `biotype`. In the [Ensembl files](http://ftp.ensembl.org/pub/current_gff3/) the `biotype`
is useful because it tells us what kind of genes they are.

Hopefully you've also downloaded GFF files for some different species and run them through your
[sqlite conversion program](Converting_gff_to_sqlite.md) to get them into a database. And, to make
this work well, hopefully you have also got your code to add an additional column taht records
which species the record came from. (I've called this column `analysis` in my code and I'll use
that below).

(For reference my version of the code is at
[solutions/part2/gff_to_sqlite.py](solutions/part2/gff_to_sqlite.py) - feel free to run that
instead, if needed.)

With this sqlite file in hand we can now start to count genes:

```
sqlite> .mode column
sqlite> .header on
sqlite> .width 50 25 25
sqlite> SELECT analysis, biotype, COUNT(*) FROM gff_data WHERE type=='gene' GROUP BY analysis, biotype ;
```
With the current data I have this gives:

    analysis                                            biotype                    COUNT(*)                 
    --------------------------------------------------  -------------------------  -------------------------
    Mus_musculus.GRCm39.104                            IG_C_gene                  13                       
    Mus_musculus.GRCm39.104                            IG_D_gene                  19                       
    Mus_musculus.GRCm39.104                            IG_J_gene                  14                       
    Mus_musculus.GRCm39.104                            IG_LV_gene                 4                        
    Mus_musculus.GRCm39.104                            IG_V_gene                  218                      
    Mus_musculus.GRCm39.104                            TEC                        3238                     
    Mus_musculus.GRCm39.104                            TR_C_gene                  8                        
    Mus_musculus.GRCm39.104                            TR_D_gene                  4                        
    Mus_musculus.GRCm39.104                            TR_J_gene                  70                       
    Mus_musculus.GRCm39.104                            TR_V_gene                  144                      
    Mus_musculus.GRCm39.104                            polymorphic_pseudogene     89                       
    Mus_musculus.GRCm39.104                            protein_coding             21834                    
    Acanthochromis_polyacanthus.ASM210954v1.104         IG_J_gene                  2                        
    Acanthochromis_polyacanthus.ASM210954v1.104         IG_V_gene                  4                        
    Acanthochromis_polyacanthus.ASM210954v1.104         TR_J_gene                  5                        
    Acanthochromis_polyacanthus.ASM210954v1.104         protein_coding             24016                    
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
    PlasmoDB-54_Pfalciparum3D7                          protein_coding_gene        5318          

This suggests [house mice](https://en.wikipedia.org/wiki/House_mouse) have about 10% more
(annotated) protein-coding genes than humans, and [spiny
chromis](https://en.wikipedia.org/wiki/Spiny_chromis) have about 10% more again.

*Note:* The various biotypes used by Ensembl are [documented
here](https://m.ensembl.org/info/genome/genebuild/biotypes.html). A `polymorphic pseudogene` is a
gene that is coding in some individuals, but not in others (including in the reference sequence).
The `IG_` and `TR_` categories are also interesting: they are the 'constant', 'joining', and
'variable' gene segments of immunuglobulin and T cell receptors. They do encode proteins, but via a
yet more [complex process that involves somatic recombination to assemble the mature genes in B and
T cells](https://en.wikipedia.org/wiki/V(D)J_recombination). These gene segments also lie in
regions that are [especially complex](https://doi.org/10.1371/journal.pcbi.1009254). However, the
vast majority of these genes are listed `protein_coding` and we will focus on these in this
tutorial.

*Note:* The above query doesn't work with the [*P.falciparum* data from
PlasmoDB](https://plasmodb.org/plasmo/app/downloads/Current_Release/)? To fix this, I kludged it by running this sql:
```
UPDATE genes SET type = 'gene', biotype = 'protein_coding_gene' WHERE analysis == 'PlasmoDB-54_Pfalciparum3D7' AND type == 'protein_coding_gene' ;
```
I wouldn't advise this type of manually-fix-the-data approach in general though, because it will break if you re-import data.

### What are all those species anyway?

If like me you're a bit unclear about what [all those
species](http://ftp.ensembl.org/pub/current_gff3/) are, now might be a good time to go and look at
the [OneZoom Tree of Life explorer](http://www.onezoom.org). This will tell you, for example, that
*Acanthochromis polyacanthus* is a member of the [Sticky Eggs
Group](http://www.onezoom.org/life/@Ovalentaria=5553750?img=best_any&anim=flight#x1307,y908,w1.5714)
 that also includes Cichlids, Silversides and Guppies, or that *mus musculus* is one of 37 species
collectively known as 'house mice'.

### Grouping data in python

The code above showed a grouping operation using SQL code. It works well, and one of its key
advantages is that it doesn't use much memory.

We could also do this in python. Indeed pandas has similar filtering and *group by* operations (see
[the pandas page on grouping](https://pandas.pydata.org/docs/user_guide/groupby.html)), just like
the above SQL code, so it is pretty easy:

```
import gff
data = gff.parse_gff3_to_dataframe("Homo_sapiens.GRCh38.104.chr.gff3" ) 
data[ data['type'] == 'gene' ].groupby( 'biotype' ).size()
```

which prints:
```
    biotype
    IG_C_gene                    14
    IG_D_gene                    37
    IG_J_gene                    18
    IG_V_gene                   145
    TEC                        1056
    TR_C_gene                     6
    TR_D_gene                     4
    TR_J_gene                    79
    TR_V_gene                   106
    polymorphic_pseudogene       49
    protein_coding            19937
    dtype: int64
````

Unfortunately however at this point we are starting to run into a potential problem: **we are using
a lot of memory**.

Whether or not you have actually run out of memory will depend on exactly how you've written your
code, and also on how much memory your machine has. Mine has 16Gb and the above uses about 4Gb of
it. This turns out to be because my [modified version of gff.py](solutions/part2/gff.py) processes
the attributes column in a way that (although it seemed a good idea at the time) uses masses of
memory.

This problem is not unusual. Because genomics data is so large, it's very easy to write code that
seems sensible and works on test data, but turns out to be a memory hog when used with real data.

[The next section](Memory_issues_and_how_to_solve_them.md) discusses this problem in more detail
and suggests ways to fix it. For now let's solve this by the simple step of not loading so much
data. I'm going to assume you have successfully created your sqlite file with some data from
different species. If so you can load just the data you need like this:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM gff_data WHERE type == 'gene'", db )
```

This dataset only has a few tens or hundreds of thousands of rows, and uses a fraction of the
memory of the full dataset. (Of course it depends on getting the data into the sqlite file first.
If memory is preventing you from doing that, I suggest trying [my version of the
code](solutions/part2/gff_to_sqlite.py), although that is still a memory hog, or if that doesn't
work, try working with smaller genomes for now). 

Now we can count genes across species in python:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM gff_data WHERE type IN ( 'gene' )", db )
data.groupby( [ 'analysis', 'biotype' ] ).size()
```

Since we're not running out of memory now, I felt happy throwing a few more species in there. This
produces:

    analysis                                      biotype               
    Acanthochromis_polyacanthus.ASM210954v1.104   IG_J_gene                     2
                                                  IG_V_gene                     4
                                                  TR_J_gene                     5
                                                  protein_coding            24016
    Camelus_dromedarius.CamDro2.104.chr.gff3      IG_C_gene                     1
                                                  IG_V_gene                    13
                                                  TR_C_gene                     1
                                                  TR_J_gene                     5
                                                  TR_V_gene                     3
                                                  protein_coding            18896
    Gallus_gallus.GRCg6a.104                      IG_V_gene                    98
                                                  protein_coding            16568
    Homo_sapiens.GRCh38.104                       IG_C_gene                    14
                                                  IG_D_gene                    37
                                                  IG_J_gene                    18
                                                  IG_V_gene                   145
                                                  TEC                        1056
                                                  TR_C_gene                     6
                                                  TR_D_gene                     4
                                                  TR_J_gene                    79
                                                  TR_V_gene                   106
                                                  polymorphic_pseudogene       49
                                                  protein_coding            19937
    Mus_musculus.GRCm39.104                       IG_C_gene                    13
                                                  IG_D_gene                    19
                                                  IG_J_gene                    14
                                                  IG_LV_gene                    4
                                                  IG_V_gene                   218
                                                  TEC                        3238
                                                  TR_C_gene                     8
                                                  TR_D_gene                     4
                                                  TR_J_gene                    70
                                                  TR_V_gene                   144
                                                  polymorphic_pseudogene       89
                                                  protein_coding            21834
    Pan_troglodytes.Pan_tro_3.0.104.chr           IG_C_gene                    11
                                                  IG_V_gene                    91
                                                  TR_C_gene                     9
                                                  TR_V_gene                    66
                                                  protein_coding            21879
    dtype: int64

Hmm... [Red junglefowl](https://en.wikipedia.org/wiki/Red_junglefowl) and [Dromedary
Camels](https://en.wikipedia.org/wiki/Dromedary) have respectively 17% and 5% fewer (annotated)
protein-coding genes than humans.

**Note.** One difference between the SQL query and the python/pandas operation, is that the pandas
one doesn't include rows with missing `biotype`. So even if you include the *P.falciparum* data
from PlasmoDB in the datat (by including fields with `type=="protein_coding_gene"`), it still won't
show in the above because the file does not record `biotypes`.

[The section on memory issues](Memory_issues_and_how_to_solve_them.md) goes into more detail about this.

## How complicated are genes?

Can we count how many transcripts each gene has?  How many exons?

To do this requires us to join the transcripts to the genes and the exons to the transcripts
somehow.  There are lots of ways of doing this - can you work out how?

**Challenge.** Write python code that reads in the appropriate data from `gff_data` table, and for
each gene counts i. the number of transcripts and ii. the average number of exons (averaged over
the transcripts for that gene).

**Hints.**

- One way to do this is to mimic the SQL approach I'll show below. If you have a table of genes, a
  table of transcripts, and a table of exons, then the [pandas merge and join
  operations](https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html#) can be used to
  link them together.
  
- Another way is to iterate through the data (for example using `.apply()`), and use it to build a data
  structure mapping genes to transcripts and transcripts to exons. For example, you might end up with something like this:
```
my_gene_summary = {
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

- Another idea is to save memory in the above by only recording the (integer) indexes of the genes and transcripts,
  instead of their (string) identifiers. How do you do that?

In the following I'll show how this can be done directly in the sqlite database and we'll then explore a pandas-based
python version.

### A sqlite approach

Just for completeness, here's how you could solve that in the sqlite database itself.

First, it's convenient to make some views of the data that are like the genes, transcripts and
exons above:

```
CREATE VIEW gene_view AS SELECT * FROM gff_data WHERE type == 'gene' ;
CREATE VIEW transcript_view AS SELECT * FROM gff_data WHERE type == 'mRNA' ;
CREATE VIEW exon_view AS SELECT * FROM gff_data WHERE type == 'exon' ;
```

Second, to make lookups efficient, we'll also need to make sure the `ID` and `Parent` fields are
indexed (and let's do `Name` as well):

```
CREATE INDEX gff_data_ID_index ON gff_data( ID ) ;
```
(This can take a minute or so as it builds the indices).

Now we can join them to make the counts.  First we'll count exons in transcripts:
```
CREATE VIEW transcript_summary_view AS
SELECT T.*, COUNT( * ) AS number_of_exons
FROM transcript_view T
LEFT JOIN exon_view E ON E.Parent == T.ID
GROUP BY T.ID;
```

If you try this out:
```
SELECT * FROM transcript_summary_view LIMIT 10 ;
```
you'll see it takes a little while to run its computation.

Next we'll summarise transcripts in genes:
```
CREATE VIEW gene_summary_view AS
SELECT G.*, COUNT(*) AS number_of_transcripts, (0.0+SUM(number_of_exons))/COUNT(*) AS average_number_of_exons
FROM gene_view G
LEFT JOIN transcript_summary_view T ON T.Parent == G.ID
GROUP BY G.ID;
```

```
SELECT * FROM gene_summary_view LIMIT 10 ;
```

We can now query it:
```
SELECT ID, Name, biotype, number_of_transcripts, average_number_of_exons  FROM gene_summary_view LIMIT 10;
SELECT * FROM gene_summary_view WHERE Name == 'ABO' ;
```

**Note.** Using `.mode line` or `.mode column` and `.header on` is useful in sqlite3. You can also
just run these queries and load them into python as described above.

**Note.** These queries are taking quite a while. In fact, we are asking sqlite to do quite a bit
of work here. This is largely because we've built this on dynamic views of the original `gff_data`
table - which is very big and includes plenty that is not relevant to our query.

One way to speed up our queries (that mimicks what we would do in python) is to 'bake' the relevant
view into their own tables first:

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

Let's use this to find the highest transcript count in each species:
```
SELECT analysis, MAX(number_of_transcripts), MAX( average_number_of_exons ) FROM gene_summary_view GROUP BY analysis ;
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

