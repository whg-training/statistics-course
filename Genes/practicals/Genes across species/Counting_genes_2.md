[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_1.md)

## How many protein-coding genes are there?

If you've followed so far you will have code `gff.py` that can parse a GFF file, and will pull out certain fields from
the `attributes` column. This includes the `ID` attribute, the `Parent` attribute that says how records are linked, and
also the `Name` column (gene names) and `biotype`. In the [Ensembl files](http://ftp.ensembl.org/pub/current_gff3/))
the `biotype` is useful and in particular, it tells us what kind of genes they are.

Hopefully you've also downloaded GFF files for some different species (if not please do this now). And, if you liked
the [sqlite part](Converting_gff_to_sqlite.md), maybe you have run them through your program `gff_to_sqlite.py` to make
a database of results.

For reference, my version of that code is in the [`solutions/part2`](solutions/part2/) folder.  Now we can do:
```
sqlite> .mode column
sqlite> .header on
sqlite> .width 50 25 25
sqlite> SELECT analysis, biotype, COUNT(*) FROM genes WHERE type=='gene' GROUP BY analysis, biotype ;
```
which gives

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

In other words [house mice](https://en.wikipedia.org/wiki/House_mouse) have about 10% more (annotated) protein-coding
genes than humans, and [spiny chromis](https://en.wikipedia.org/wiki/Spiny_chromis) have about 10% more again. 

*Note:* The various biotypes used by Ensembl are [documented
here](https://m.ensembl.org/info/genome/genebuild/biotypes.html). The `IG_` and TR_` categories are
interesting: they are the 'constant', 'joining', and 'variable' gene segments of immunuglobulin and
T cell receptors. They do encode proteins, but via a yet more [complex process that involves
somatic recombination to assemble the mature genes in B and T
cells](https://en.wikipedia.org/wiki/V(D)J_recombination). These gene segments also lie in regions
that are [especially complex](https://doi.org/10.1371/journal.pcbi.1009254). (But we will focus on
regular protein-coding genes in this tutorual.)

*Note:* Does the above query work with the [*P.falciparum* data from PlasmoDB]()? How many
protein-coding genes does *P.falciparum* have?

### Grouping data in python

The code above showed a grouping operation using SQL code. It works well. One of its key advantages
is that it doesn't use much memory.

We could also do this in python. In fact pandas has similar filtering and *group by* operations
(see [the pandas page on grouping](https://pandas.pydata.org/docs/user_guide/groupby.html)), just
like the above SQL code, so it is pretty easy:

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

Unfortunately however at this point we are starting to run into a possible problem: **we are using
a lot of memory**.

Whether or not you actually run out of memory will depend on exactly how you've written your code,
and also on how much memory your machine has. Mine has 16Gb and the above uses about 4Gb of it.
This turns out to be because my [modified version of gff.py](solutions/part2/gff.py) processes the
attributes column in a way that (although it seemed a good idea at the time) uses masses of memory.

This problem is not unusual. Because genomics data is so large, it's very easy to write code that
seems sensible and works on test data, but turns out to be a memory hog when used with real data.

[The next section](Memory_issues_and_how_to_solve_them.md) dicsusses this problem in more detail
and suggests ways to fix it. For now let's solve this by the simple step of not loading so much
data. I'm going to assume you have successfully created your sqlite file with some data from
different species. If so you can load just the data you need like this:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM genes WHERE type == 'gene'", db )
```

This dataset only has a few tens or hundreds of thousands of rows, and uses a fraction of the memory of the full
dataset. (Of course it depends on getting the data into the sqlite file first. (If memory is preventing you from doing
this, one option is to work with smaller genomes for now). We can now count genes across species in python:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM genes WHERE type IN ( 'gene' )", db )
data.groupby( [ 'analysis', 'biotype' ] ).size()
```
