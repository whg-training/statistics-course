[Up to table of contents](Introduction.md)
[Back to the previous page](Converting_gff_to_sqlite.md)

## Gathering some first statistics

Let's start to attack our scientific questions. From simplicty, from here on in we will work with the GFF3 files
from the [Ensembl ftp site](http://ftp.ensembl.org/pub/current_gff3/) - also those from
[PlasmoDB](https://plasmodb.org/plasmo/app/downloads/Current_Release/). These use the terminology `mRNA` for a gene
transcript, and they also have the genome sequence lengths written in the metadata, making life easy for us.

We will focus on protein-coding genes, and their transcripts, exons and coding sequence. They have type `gene`, `mRNA`, `exon`
and `CDS` in the files respectively. The basic setup is:

- Each gene can have multiple transcripts (i.e. multiple expressed forms at mRNA level - e.g. they might differ in how
  exons are spliced together, or have different transcription start or end sites).
  
- Each transcript is made up of one or more exons. The introns in between are spliced out. (This happens during
  transcription process.)

- Of the exon sequence that makes it into a transcript, only a portion is actually translated into the mature protein
  during processing by the ribosomes.  This is the *coding sequence* and the rest is *untranslated sequence*.

What's in these files is humanities' best guess at what this picture looks like in each organism.

### 

We're actually in a good shape for many of these questions.  Let's take them one by one:

### How many genes are there?

Let's start with counting genes - easy!  For example using the sqlite file you created in the last step:

```
sqlite> SELECT COUNT(*) FROM genes WHERE type==`gene` ;
```

(I ran this interactively in the sqlite3 shell.  Run `sqlite3 genes.sqlite` to get there.  Type `.quit` to quit.)

If you used my code or followed the suggestions, it adds an `analysis` column so you could do this for multiple species:

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

The above counts all "genes" - not all of these code for proteins though.  

**Question:** How do we count protein-coding genes?

**Challenge:** Can you update your `gff.py` code to extract other useful fields from the `attributes` column? In
particular - the `biotype` field is relevant. (And why don't we also get the `Name` field, so we can easily see the
gene name).

My version of that is in the [`solutions/part2`](solutions/part2/) folder.  Now we can do:
```
sqlite> .mode column
sqlite> .header on
sqlite> .width 50 25 25
sqlite> SELECT analysis, biotype, COUNT(*) FROM genes WHERE type=='gene' GROUP BY analysis, biotype
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

*Question:* how do you do this type of counting in python / pandas?  For one GFF file you could do something like this:
```
import gff.py
data = gff.parse_gff3_to_dataframe("Homo_sapiens.GRCh38.104.chr.gff3" ) 
data.
```

In other words [House mice](https://en.wikipedia.org/wiki/House_mouse) have about 10% more (annotated) protein-coding
genes than [Humans](), and [Spiny chromis](https://en.wikipedia.org/wiki/Spiny_chromis) have about 10% more again. 

*Task:* Run this on more species - either by converting to sqlite (*Warning*: keep an eye on the file size) or
loading the file directly using `gff.py`.

*Note:* The various biotypes [used by Ensembl are [documented
here](https://m.ensembl.org/info/genome/genebuild/biotypes.html). The `IG_` and TR_` categories are particularly
interesting because they encode the 'constant', 'joining', and 'variable' gene segments of immunuglobulin and T cell
receptors. They do encode proteins, but via a yet more complex process that involves somatic recombination to assemble
the mature genes in B and T cells. If you want to learn about the complexity of these regions, check out [this recent
paper by Jia-Yuan Zhang](https://doi.org/10.1371/journal.pcbi.1009254). (But we will focus on protein-coding genes here.)

*Note:* Does the above query work with the [*P.falciparum* data from PlasmoDB]()?  How many protein-coding genes does *P.falciparum* have?


