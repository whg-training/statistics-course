## A basic python approach

We could certainly start by reading / parsing each line
of the file, ignoring the complicated relationship-between-record-y bits. That should be easy right? Something like this:

```
def readGFF3( filename ):
    """Read a GFF3-formatted file, output (something)"""
    # do something
    pass
````

**Note: ** If you are not familiar with python syntax, now would be a good time to refresh via any of the available tutorials.
The above code defines a function, and currently contains a documentation comment (the `"""..."""` bit) and also a code
comment (starting with `#`.). The `pass` is just there to make this a valid python function that does nothing at the
moment.  I'll paste equivalents in other languages below.

Can you write this? 

But hang on... *Should* you write this? As the comment indicates, there're already a couple of issues to deal with:

**Issue 1.** What should `readGFF3()` take in as a parameter?

The best way to figure this out is to think about how you will call it from other code. The seemingly obvious choice above is
to make it take in a filename, so we would call it like this:

```
readGFF3( "/path/to/gff/file.gff" )
```

That'd be ok. But hold on. Our principles say we should write easily testable code. To test the function above, we'd
have to first put test data into a file on the filesystem, then pass the path in... too complex! This suggests we
should write our function to take in some data instead - let's say, an array of lines of data. So then we can test it
easily:

```
testData = [
    "##gff-version 3",
    |#description: test data",
    "chr1	me	gene	1	1000	.	+	.	ID=my_test_gene;gene_id=my_test_gene.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2",
    "chr1	me	exon	1	1000	.	+	.	ID=my_test_exon;gene_id=my_test_exon"
]
result = readGFF3( testData )
# do something to check the result here
```

But hang on - if we're not actually *reading* the data (from a file) then the function name is wrong. It is only parsing
the data. So we should actually be writing this function instead:
```
def parseGFF3Data( data ):
    """read an array of lines of GFF3-formatted data, output (something)"
    ...
```

**Conclusion:** we shouldn't write `readGFF3()`.  We should write `parseGFF3Data()` instead.

**Consideration 2.** What should `parseGFF3Data()` output?

There are basically two 

The 'keep it simple' principle pretty much dictates how we do this:

- The GFF file is made up of multiple data rows - this says we should return an array (that is, a python `list`) of
  rows.
  
- The [GFF spec](https://m.ensembl.org/info/website/upload/gff3.html) makes it pretty clear that each row has 9
  columns.
  
- The first 8 are single values that contain either a string, an integer (`start`, `stop`), or a floating-point number
  (`score`). These map to python datatypes.

- The last column (`attributes`) is a kind of catch-all that contains key-value pairs. The
  simplest way to represent a set of key-value pairs in python is a `dict`.

This suggests the following structure:

```
class GFFRecord:
    def __init__( self ): # This is a python 'constructor'
        self.seqid = None
        self.source = None
        self.type = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = none
        self.phase = None
        self.attributes = {}  # a python dict
```

And now we can write our function:
```
def parseGFF3Data( data ):
    """Given an array of lines of GFF3-formatted data,
    parses the data and outputs an array of GFFRecord structs."""
    result = []
    # fill in result here
    return result
```

*Note.* Arguably, an even simpler way is to forget the above class and just return a dict (or even an array, that
wouldn't have any column names) for each row. I've gone with the above because I think it makes *using* the rows
simpler (as opposed to parsing them).

### Thoughts ###

We haven't written any code that does anything yet. But a bit of reasoning has led us to code that is better
than the code we would have written (or at least the code *I* would have written) in a couple of ways:

- `parseGFF3Data()` will be shorter and simpler to write than `readGFF3()` - since it doesn't have to worry about opening files etc.
- `parseGFF3Data()` is easier to test than `readGFF3()`, since we can pass data directly in.
- `parseGFF3Data()` is arguably easier to understand (because it is better named.
- It's clear what `parseGFF3Data()` is supposed to return.

**Note.** It's never that easy in practice, and you probably won't reach the nicely-testable-well-written code first
time in a real problem. However, you can get there by rewriting (or 'refactoring') your code each time you visit it
with an eye principles like the ones above.

We can test it now using our test data above:

```
testData = [
    "##gff-version 3",
    "#description: test data",
    "chr1	me	gene	1	1000	.	+	.	ID=my_test_gene;gene_id=my_test_gene.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2",
    "chr1	me	exon	1	1000	.	+	.	ID=my_test_exon;gene_id=my_test_exon"
]

result = parseGFF3Data( testData )
assert result.length == 2
assert result[0].seqid == 'chr1'
assert result[0].type == 'gene'
assert result[0].attributes['ID'] == 'my_test_gene'
assert result[1].attributes['ID'] == 'my_test_exon'
# ... etc.
```

Go test-driven development!  You can run this code right now, but of course it fails because we haven't written the actual code.

## Writing `parseGFF3Data()`

If you follow the same appraoch as above, it's now pretty easy to write `parseGFF3Data()`. The data comes in rows, so
we should iterate over them and parse them one by one:

def parseGFF3Data( data ):
    """Given an array of lines of GFF3-formatted data,
    parses the data and outputs an array of GFFRecord structs."""
    result = []
    for line in data:
        result.append( parseGFF3Record( line ))
    return result

```
Solved!



If you get here

To make it work, you need to to [write the parseGFF3Data function][Writing_the_parseGFF3Data_function.md].
