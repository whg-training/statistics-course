[Up to table of contents](Introduction.md)
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

*Note:* The various biotypes [used by Ensembl are [documented
here](https://m.ensembl.org/info/genome/genebuild/biotypes.html). The `IG_` and TR_` categories are particularly
interesting because they encode the 'constant', 'joining', and 'variable' gene segments of immunuglobulin and T cell
receptors. They do encode proteins, but via a yet more complex process that involves somatic recombination to assemble
the mature genes in B and T cells. If you want to learn about the complexity of these regions, check out [this recent
paper by Jia-Yuan Zhang](https://doi.org/10.1371/journal.pcbi.1009254). (But we will focus on protein-coding genes here.)

*Note:* Does the above query work with the [*P.falciparum* data from PlasmoDB]()?  How many protein-coding genes does *P.falciparum* have?

### Grouping data in python

The code above showed a grouping operation using SQL code. It works well and one of its key advantages is that
it doesn't use much memory.

We could also do this in python. In fact pandas has filtering and *group by* operations (see [the pandas page on
grouping](https://pandas.pydata.org/docs/user_guide/groupby.html)), just like the above SQL code, so it is pretty easy:

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

(**Note:**  you might need to first run `pandas.set_option('display.max_rows', 20 )` to see all the rows above.)

### Memory issues

Unfortunately we are starting to run into a major problem: **we are running out of memory**.

When I run the above three lines with the [rewritten gff.py](solutions/part2/gff.py) in a freshly-started python instance,
the process uses lots of memory - about 4Gb.  And my computer has only 16Gb in total!  If you were to try to
analyse large data - say multiple species at once - you'd quickly find all your memory used up.

This type of problem is actually quite typical for bioinformatics analyses - they always get larger until we hit some
limit. High-level approaches like the one we've taken (which loads all the data into memory and then processes it) seem
good at first, but they do not control memory usage.  And that rapidly becomes a problem as data volumes grow.

If you have hit this problem then it is a pain! In general the main ways to solve it are:

- work with data subsets (genome regions, specific record types, smaller organisms).
- use memory-efficient data structures (like the sqlite database we are using).
- use things like sqlite that allow efficient access to on-disk data.
- write code more carefully to control memory usage.

[The next section](Memory_issues_and_how_to_solve_them.md) goes into more detail about this.

What's more, with these data volumes, seemingly innocuous changes to this type of code can start to make a huge
difference. My two versions of `gff.py` - [version 1](solutions/part1/gff.py) and [version 2](solutions/part2/gff.py) -
do almost the same thing. The second one just adds a couple of columns. But it turns out to use about 4 times the
memory.

Not surprisingly, in our program this is all to do with attribute parsing. They key difference between the two versions
is that the original version of `parse_gff3_to_dataframe()` **never stored the parsed/unpacked attributes strings**.
(Instead - [as profiling revealed](Converting_gff_to_sqlite.md) - it was wasting time by parsing them twice). But the
second version does do this: on line 25 it says:

```
    attributes = result['attributes'].apply( parse_attributes )
```
This one line turns out to use up about 2.5Gb of memory!

Considerations like this mean we need to add an additional aim when we are coding:

- it ought to work
- it ought to not take too long to do it
- it ought to be obvious what it does
- **it ought not to waste memory**

**Question.** How much memory is your version of the code using?

The main approaches to solve this problem are:

- work with data subsets (genome regions, specific record types, smaller organisms)
- use memory-efficient data structures (like the sqlite database we are using)
- write code more carefully to control memory usage

**Advanced challenge:** write a version of `gff_to_sqlite.py` that has low memory usage.

**Hints**:

- There's really no need for `gff_to_sqlite.py` to use any memory at all! It could process rows one at a time instead
  of loading them all in.

- The [the pure python version](gff_to_sqlite_python_version.py) is a good place to look for the needed sql statements.

- To make this slick, you should aim to commit the rows to the database in chunks (say of 10,000 rows). Then call
  `db.commit()` to write the data to the file after every chunk.

- Don't forget to commit the last chunk.

## Grouping data in python - revisited

As an example of working with subsets - let's just load the subset of records we need:

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
