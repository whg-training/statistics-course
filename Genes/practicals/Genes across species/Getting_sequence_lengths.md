[Up to table of contents](README.md)
[Back to the previous page](Memory_issues_and_how_to_solve_them.md)

## Getting sequence lengths

To answer questions like 'how much of the genome is in genes' we need to know what the length of each genome is.
Happily, if you look at the [GFF files on Ensembl](http://ftp.ensembl.org/pub/current_gff3/) the sequences from the
relevant reference are listed in the metadata. For example, for humans they look like this:

    ##gff-version 3
    ##sequence-region   1 1 248956422
    ##sequence-region   10 1 133797422
    ##sequence-region   11 1 135086622
    ##sequence-region   12 1 133275309
    ...
    ##sequence-region   7 1 159345973
    ##sequence-region   8 1 145138636
    ##sequence-region   9 1 138394717
    ##sequence-region   MT 1 16569
    ##sequence-region   X 1 156040895
    ##sequence-region   Y 2752083 56887902

The three pieces of data on each row are: the sequence name (first column), the sequence starting position, and the
ending position. (Incidentally, the reason the Y chromosome starts at a base different than 1 is because it omits the
[pseudoautosomal regions](https://www.ncbi.nlm.nih.gov/grc/human?filters=chr:Y#current-regions), which are often
treated for analysis purposes as being on the X chromosome).

If you feel like it, please try the following:

**Challenge:** write a function `parse_sequences_from_gff_metadata()` that

- takes the gff filename as argument
- return a pandas dataframe with the above information.  It should have three columns, the `sequence`, the `start`, and the `end`.

Use this function in `gff_to_sqlite.py` so that it records sequence lengths into the output
database as well.

Here's a test case:

```
import unittest

class TestGff(unittest.TestCase):
    def test_parse_sequences_from_gff_metadata():
       import io
       test_data = """##gff-version 3
##sequence-region   1 1 248956422
##sequence-region   10 1 133797422
1       GRCh38  chromosome      1       248956422       .       .       .       ID=chromosome:1;Alias=CM000663.2,chr1,NC_000001.11
"""

       sequences = parse_sequences_from_gff_metadata( io.StringIO( test_data ))

       self.assertEqual( sequences['sequence'][0], '1' )
       self.assertEqual( sequences['sequence'][0], '10' )
       self.assertEqual( sequences['start'][1], 1 )
       self.assertEqual( sequences['start'][0], 1 )
       self.assertEqual( sequences['end'][0], 248956422 )
       self.assertEqual( sequences['end'][1], 133797422 )
```

If you want some hints or just don't fancy writing this piece of sequence parsing code - feel free
to look at [my version is here](solutions/part3/gff.py). And I added it to
[`gff_to_sqlite.py`](solutions/part3/gff_to_sqlite.py).

We are now ready to [figure out how much of the genome is in genes](How_much_of_the_genome_is_in_genes.md).
