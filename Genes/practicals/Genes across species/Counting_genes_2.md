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


