[Up to table of contents](README.md)

[Back to the previous page](Getting_sequence_lengths.md)

[Go to the next page](Visualiation.md)

## How much of the genome is in genes?

To figure out how much of the genome is covered by genes, or by exons, we face a problem.
In principle we could just add together the gene lengths.
```
genes['length'] = genes['end'] - genes['start']
genes.groupby('analysis').agg(
    covered_length = pandas.NamedAgg( column = 'length', aggfunc = lambda w: w.sum() )
)
```

Then we could compare to the sequence lengths
```
sequence_lengths = pandas.read_sql( "SELECT analysis, SUM(end-start) AS sequence_length FROM sequences GROUP BY analysis", db )
```

So this would say for example that about 44% of the human genome is covered by genes.

### Union of regions

Unfortunately this isn't so simple - because many genes overlap each other. This happens either
because there genuinely are different genes encoded by the same bit of DNA, or because of
additional annotated 'genes' that arise due to the computational gene annotation process. (We saw
one earlier - [a tiny annotated gene inside *SLC4A2*](Counting_genes_3.md).)

To handle this we have to compute the regions covered by a bunch of overlapping genes, and for that
we need to be able to compute this overlap.

**Challenge.**. Write a function `compute_union_of_regions()` that computes the total region covered
by a set of (possibly overlapping ranges). Both the input and output should be a list of pairs of
the form `[ start, end ]` (where `end` >= `start` and the coordinates are all non-negative
integers. The output should contain a set of non-overlapping ranges that cover the same set of
positions as the input ranges. And for testability reasons, let's also require the output to be
sorted (by the region start position).

Here is a test:
```
class TestRanges(unittest.TestCase):
    def test_union_of_ranges( self ):
        # check some simple cases first
        self.assertEqual( compute_union_of_regions( [[1,10]] ) == [[1,10]] )
        self.assertEqual( compute_union_of_regions( [] ) == [] )
        self.assertEqual( compute_union_of_regions( [[1,10], [11,11]] ) == [[1,11]] )

        test_data = [
            [1, 10],
            [19,199],
            [5, 6],
            [9, 15],
            [20, 25]
        ]
        result = compute_union_of_ranges( test_data )
        self.assertEqual( len( result ), 2 )
        self.assertEqual( result[0] = [ 1, 15 ] )
        self.assertEqual( result[1] = [ 19, 15 ] )

```

**Hints.** 

* First [sort the list of regions by the start
  point](https://docs.python.org/3/howto/sorting.html). (But be aware that python functions can
  mutate their arguments). You may want to use [the `sorted()` function]() rather than sorting
  in-place.
  
* Now traverse the list of regions, keeping track of the current interval and extending it if
  necessary when you encounter overlapping input regions.

**Important note.** The coordinates in the GTF file are defined to [follow the 1-based
convention](http://www.ensembl.org/info/website/upload/gff.html). This means that the genome
coordinates start at 1, and also that regions are defined to be closed - i.e. they contain both
their endpoints. A region like [1,10] therefore contains 10 base positions.

(If this sounds obvious, I'm raising it because in other contexts a 0-based, half-open convention is
used instead (in which the region [1,10) would contain only 9 positions, and would miss the 1st
genome location at zero). This is true for [the database that underlies the [UCSC Genome
Browser](https://genome-blog.soe.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
 for example (although not the browser itself, which converts coordinates to 1-based), and is
common in programming generally.)

My solution is [here](solutions/part3/gff/regions.py).

### So how much of the genome is in genes?

We're ready to answer this! Let's use the above to write function that counts the bases covered in each genome:

```
# Add this to your gff.py file
def compute_genome_bases_covered( regions, sequences ):
   """Given a set of regions (as a dataframe with analysis, seqid, and start, and end columns),
   and a set of sequences (as a dataframe with analysis, seqid and sequence_length columns, return
   a dataframe showing the total number and proportion of sequence bases covered by the regions in
   each analysis."""
   def sum_region_lengths( regions ):
      # We convert from pandas to a list...
      aslist = regions[['start', 'end']].values.tolist()
      # ...compute the union using your function...
      union = compute_union_of_regions( aslist )
      # ...and sum to get the result.
      result = sum( [ (region[1]-region[0]+1) for region in union ] )
      # Finally, we return as a pandas.Series() object
      # because this lets us name the output variable.
      return pandas.Series(
         sum( [ (region[1]-region[0]+1) for region in union ] ),
         index = ['bases_covered']
      )
   def compute_bases_covered_in_one_chromosome( regions ):
      return (
         regions
         .groupby( [ "analysis", "seqid" ])
         .apply( sum_region_lengths )
      )
   def add_sequence_lengths( coverage, sequences ):
      return pandas.merge(
         coverage,
         sequences,
         left_on = [ "analysis", "seqid" ],
         right_on = [ "analysis", "seqid" ]
      )
   def sum_over_chromosomes( coverage ):
      result = coverage.groupby( "analysis" ).agg(
         bases_covered = pandas.NamedAgg( column = "bases_covered", aggfunc = sum ),
         sequence_length = pandas.NamedAgg( column = "sequence_length", aggfunc = sum )
      )
      result['proportion'] = result['bases_covered'] / result['sequence_length']
      return result
   per_chromosome = add_sequence_lengths(
      compute_bases_covered_in_one_chromosome( regions )
   )
   return sum_over_chromosomes( per_chromosome )
```

The code above is a bit complex, but in short:

* it uses the `compute_union_of_regions()` function to compute the sum of lengths of genes *in each
  chromosome and species*.
* it then adds a column with the sequence lengths [that we computed earlier](Getting_sequence_lengths.md) to the result.
* Finally it summarises over the whole genome for each species and returns the result.

Let's try it out:
```
import pandas, sqlite3, gff
db = sqlite3.connect( "genes.sqlite" )
genes = pandas.read_sql( """
   SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype
   FROM gff_data
   WHERE type == 'gene'
   AND biotype == 'protein_coding'
""", db
)
sequences = pandas.read_sql( "SELECT analysis, seqid, ( end - start ) AS sequence_length FROM sequences", db )
gff.compute_genome_bases_covered( genes, sequences )
```

This gives:

                                                bases_covered  sequence_length  proportion
    analysis                                                                               
    Acanthochromis_polyacanthus.ASM210954v1.104      415285145        830196730    0.500225
    Camelus_dromedarius.CamDro2.104.chr.gff3         815783547       2052758671    0.397408
    Gallus_gallus.GRCg6a.104                         480522985       1050156563    0.457573
    Homo_sapiens.GRCh38.104                         1298183323       3085194781    0.420778
    Mus_musculus.GRCm39.104                         1046654226       2723431121    0.384315
    Pan_troglodytes.Pan_tro_3.0.104.chr             1101889770       2967125051    0.371366
    PlasmoDB-54_Pfalciparum3D7                        18174434         23332823    0.778921

So only 42% of the human genome is covered by protein-coding genes. This drops a bit in
chimpanzees, while a colossal 78% of the *P.falciparum* genome is in genes.

### How much of the genome is in exons?

Here we'd better be a bit careful to only get exons in protein-coding genes. Let's do it in python
as follows (as usual this could be done in sql as well, if we wanted):

```
transcripts = pandas.read_sql( """
   SELECT analysis, ID, Parent, seqid, start, end
   FROM gff_data WHERE type == 'mRNA'
""", db )

exons = pandas.read_sql( """
   SELECT analysis, ID, Parent, seqid, start, end
   FROM gff_data WHERE type == 'exon'
""", db )

print( "%d transcripts, %d exons.\n" % ( transcripts.shape[0], exons.shape[0] ))
```

Now find just the transcripts that actually are protein-coding genes:

```
transcripts = transcripts[ transcripts['Parent'].isin( genes['ID'])]
exons = exons[ exons['Parent'].isin( transcripts['ID'] ) ]
print( "%d transcripts, %d exons.\n" % ( transcripts.shape[0], exons.shape[0] ))

gff.compute_genome_bases_covered( exons, sequences )
```

This gives:

                                                 bases_covered  sequence_length  proportion
    analysis                                                                               
    Acanthochromis_polyacanthus.ASM210954v1.104       64427980        830196730    0.077606
    Camelus_dromedarius.CamDro2.104.chr.gff3          44585999       2052758671    0.021720
    Gallus_gallus.GRCg6a.104                          41635090       1050156563    0.039647
    Homo_sapiens.GRCh38.104                           85003514       3085194781    0.027552
    Mus_musculus.GRCm39.104                           76703787       2723431121    0.028164
    Pan_troglodytes.Pan_tro_3.0.104.chr               61897466       2967125051    0.020861
    PlasmoDB-54_Pfalciparum3D7                        16683959         23332823    0.715042

For most genomes, even though a large proportion of genome is in genes, hardly any is in exons. But
P.falciparum is different: we already saw it has [lots of single-exon genes](Counting_genes_3.md)
and sure enough most of the genome is in these exons.

### How much of the genome is in coding sequence?

This is easy now:

```
cds = pandas.read_sql( """
   SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype
   FROM gff_data WHERE type == 'CDS'
""", db)
cds = cds[ cds["Parent"].isin( transcripts['ID'])]
gff.compute_genome_bases_covered( cds, sequences )
```

This gives:

                                                 bases_covered  sequence_length  proportion
    analysis                                                                               
    Acanthochromis_polyacanthus.ASM210954v1.104       37182784        830196730    0.044788
    Camelus_dromedarius.CamDro2.104.chr.gff3          31455083       2052758671    0.015323
    Gallus_gallus.GRCg6a.104                          28472131       1050156563    0.027112
    Homo_sapiens.GRCh38.104                           35837191       3085194781    0.011616
    Mus_musculus.GRCm39.104                           36746474       2723431121    0.013493
    Pan_troglodytes.Pan_tro_3.0.104.chr               35836571       2967125051    0.012078
    PlasmoDB-54_Pfalciparum3D7                        12357002         23332823    0.529597

In humans, only about half of the sequence in protein-coding gene exons actually codes for proteins - genes have a great deal of *untranslated sequence*.

### Conclusions

Two really striking things jump out from this to me.  First, in many eukaryotic organisms, hardly any of the
genome is actually coding for protein. 

What is the rest doing? Well, we saw that some of it is in [pseudogenes or
encodes non-protein-coding RNAs](Counting_genes_1.md). Quite a lot actually is functional [as
described in the ENCODE paper](https://doi.org/10.1038/nature11247) (at least in the limited sense that
proteins bind there, as part of a general process of regulating gene expression). Some of it is viral DNA
that has been retro-transposed into the genome - such as the [LINEs and SINEs in the human
genome](https://en.wikipedia.org/wiki/Retrotransposon).  Some is genuinely doing nothing. However, much
of the function of even well-annotated genomes (like humans) is not known.

On the other hand, some organisms like malaria seem to be making use of a great deal of their
genome for making proteins. Quite why this is is also not known.

## Where next?

Here are some [suggested next steps](Where_next.md).
