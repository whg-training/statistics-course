[Up to table of contents](README.md)
[Back to the previous page](Getting_sequence_lengths.md)

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
sequence_lengths = pandas.read_sql( "SELECT analysis, SUM(end-start) FROM sequences GROUP BY analysis", db )
```

So this would say that about 44% of the human genome is covered by genes.

Unfortunately this isn't so simple **because many genes overlap each other**. We have to compute the regions covered by a
bunch of overlapping genes, and for that we need interval arithmetic.

**Challenge.**. Write a function `compute_union_of_ranges()` that computes the total region covered by a set of
(possibly overlapping ranges). Both the input and output should be a list of pairs of the form `[ start, end ]` (where
`end` >= `start` and the coordinates are all non-negative integers). The output should contain no overlapping ranges,
so that the total region covered can be computed by summing the region lengths.

Here is a test:
```
class TestRanges(unittest.TestCase):
    def test_union_of_ranges_( self ):
        # check some simple cases first
        self.assertEqual( compute_union_of_ranges( [[1,10]] ) = [[1,10]] )
        self.assertEqual( compute_union_of_ranges( [] ) = [] )

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

My solution is [here](..)
