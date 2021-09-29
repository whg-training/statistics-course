## Writing some code to process GFF files

Your task, if you choose to accept it, is to write some reuseable code that processes a GFF file,
make sense of the data, and gathers some statistics. And then to apply it to analyse genes from
multiple species, to hopefully learn something about the genome biology  of the world's organisms.

We would like **you** to write the code that does this - you can either do this on your own or
working together in a group if you prefer. There will be lots of support, and this tutorial will
also guide you through one way to do it. We'll then talk about both the results and the code itself
in the discussion session.

We're not writing this code for it's own sake but to answer our questions like the ones in our introduction:

- How many genes are there?
- How big are they?
- How much of the genome is in genes?
- How complex are genes - How many exons?  How many different transcripts?
- How much of genes is actually protein-coding sequence - and how much is untranslated?
- How much do these patterns differ across species?

How to write this code?  Well there are a few ways:


**Do it yourself.**. It may well be that you already have a good idea how to go about this. If so,
feel free to dive straight in. You're free to use any language or system you like for this -
standard options might be [python](https://www.python.org) or [R](https://cran.r-project.org), but
you could also use [julia](https://julialang.org), or even
[C++](https://en.wikipedia.org/wiki/C%2B%2B) or another compiled language. 

**Use a package.** You are also free to use whatever packages you like. Chances are you will find a
package that will parse this data for you. And that's fine.

**Follow this tutorial.** For top marks you have to write your own code! This tutorial will do it
in python, either relying on `pandas` data frame library, or writing it by hand.

## Diving straight in

If you [looked at the gene annotation data](What_gene_annotation_data_looks_like.md), you'll know
that it comes in rows of data that are tab-delimited, but that it is also kind of *relational*
(meaning that the records refer to each other, via the `Parent` attribute). Somehow we have to get
this into a data structure that we can work with.

According to [our principles](Introduction.md), we should keep it simple. So let's do the simplest
thing possible and start with the bit that parses the data (ignoring the hierarchical structure for
now.) That should be easy right?  Here's a python function that does it:

```
def readGFF3( data ):
    """Read GFF3 data, output (something)"""
    result = <something here>
    # do some work.
    return( result )
````

**Note:** If you are not familiar with python syntax, now would be a good time to refresh via any of the available tutorials.
The above code defines a function, and currently contains a documentation comment (the `"""..."""` bit) and also a code
comment (starting with `#`.). The `pass` is just there to make this a valid python function that does nothing at the
moment.  I'll paste equivalents in other languages below.

Right away however we must make some choices. What should `readGFF3()` return? And what argument should it take?

If you want to fully test-driven, you could now write a test that specifies this - even before we've written any code.
For example we could write this:
```
testData = """##gff-version 3
#description: test data
chr1	me	gene	1	1000	.	+	.	ID=my_test_gene;gene_id=my_test_gene.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2
chr1	me	exon	1	1000	.	+	.	ID=my_test_exon;gene_id=my_test_exon
"""
result = readGFF3( testData )
assert result[0].type == 'gene'
assert result[0].start == 1
assert result[1].attributes['ID'] == 'my_test_exon'
```
This dictates that the function had better be able to take in a string of data, and it had better
return an array of objects that have properties we can access.

On the other hand - it's often better to think of datasets in terms of data frames - i.e. as great
big tables of data, with many rows and a fixed number of columns that have particular names and
data types. The GFF3 format is like this - it's just a big table with N rows and 9 columns
(although one of them is that awkward `attributes` list of key/value pairs.) So this view would say
that `readGFF3()` should return a data frame:

```
result = readGFF3( testData )
assert result.shape == (2,11)
assert result['type'][0] == 'gene'
assert result['start']0] == 1
```

This looks pretty similar to the above but now we think of the function as returning a
2-dimensional data frame with a definite shape (i.e. 2 rows and 9 columns).




Can you write this function?

But hang on... *Should* you write this? As the comment indicates, there're already a couple of issues to deal with:

**Issue 1.** What should `readGFF3()` take in as a parameter?

The best way to figure this out is to think about how you will call it from other code. The seemingly obvious choice is
to make it take in a filename, so we would call it like this:

```
readGFF3( "/path/to/gff/file.gff" )
```

That looks sensible. But wait, our principles say we should make our code easy to test. To test the above function
above we'd have to first put test data into a file on the filesystem, and pass the path in - which is complicated. Or
else to somehow trick the function into thinking it has a file.  This is of course doable
Instead this principle suggests we should write
our function to take in some data instead - let's say, an array of lines of data. So then we can test it easily:

```
testData = """##gff-version 3
#description: test data
chr1	me	gene	1	1000	.	+	.	ID=my_test_gene;gene_id=my_test_gene.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2
chr1	me	exon	1	1000	.	+	.	ID=my_test_exon;gene_id=my_test_exon
"""
result = readGFF3( testData )

# Check it was right here...
```

Going a bit further, if we're not actually *reading* the data (from a file) in the function then the function name is wrong. It is only parsing
the data, so we should rename it:
```
def parseGFF3Data( data ):
    """read GFF3-formatted data, output (something)"""
    ...
```

**Conclusion:** we shouldn't write `readGFF3()`.  We should write `parseGFF3Data()` instead.

**Consideration 2.** What should `parseGFF3Data()` output?

There are basically two approaches we could take here. The first and probably easiest is to treat the data as a data
frame. This means we get to use all of the data frame apparatus (i.e. the `pandas` library in python) which makes life
easy. Our function is:

```
def parseGFF3Data( data ):
    """read GFF3-formatted data, output a data frame of the parsed data"""
    ...
```

If you explore [pandas](https://pandas.pydata.org/docs/) a little you'll find the `read_table` function which does most
of the work for us:

```
def readGFFFile( filename ):
    import pandas
    return pandas.read_table(
        filename,
        comment = '#',
        names = [ 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )

def parseGFF3Data( data ):
    """read GFF3-formatted data, output a data frame of the parsed data"""
    readGFFFile( data )
    return result
```

However, this isn't good enough for a couple of reasons:

- The data contains missing values that appear in the data as `.`.  They should come out as missing values (`None`) in python.
- The rows contain important `ID` and `Parent` columns that are buried in the `attributes` column.  We would like these as first-class 

 input data is basically rectangular (many rows by 9
columns). This is what a *data frame* exists to model. Therefore the function should return a data frame. This is what
most people want to do because (as long as you don't mind using the `pandas` library it gives an easy-to-use solution.)

## A lower-level appraoch

There's no particular need to use data frames or `pandas` here, and actually the GFF3 file format more or less dictates what the
function could do:

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

We could certainly start by ignoring the difficult bits and just reading / parsing the file, ignoring the complicated
relationship-between-record-y bits. That should be easy right? Something like this:

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

But hang on... *Should* you write this? As the comment indicates, there're already a couple of issues to deal with

### Loading the data

If you are a glutton for punishment you could try doing this [the low-level way](parsing_gff3_the_low_level_way.md).
However, like many bioinformatics data formats, these files are basically 2x2 tables, with a number of rows each with a
fixed number of columns (that have names).  In other words it is a **dataframe**.  So we want something like this:

```
def readGFF3( filename ):
    """Read GFF3-formatted data from the given filename, return a dataframe"
```

This is what data frames were invented for. The simplest way to load in
python is to use the `pandas` data frame library:
```
import pandas
gffData = pandas.read_table(
    "/path/to/gff/file.gff",
    comment = '#',
    names = [ 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
)
print(gffData)
```

Who wants to type that every time?  Not me, so let's wrap it into a function:
```
def readGFF3FileToDataframe( filename ):
    import pandas
    return pandas.read_table(
        filename,
        comment = '#',
        names = [ 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
```

Go ahead and try it:
```
readGFF3FileToDataframe(
    "PlasmoDB-54_Pfalciparum3D7.head.gff"
)
```

**Note.** Pandas takes a do-everything-for-you approach, which can be quite helpful. (On the other hand
arguably goes against the coding principles I was suggesting). For example, it will happily take a URL:
```
    readGFF3FileToDataframe(
        "http://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/PlasmoDB-54_Pfalciparum3D7.gff"
    )
```
When I try this I get an error to do with https and certificates, though.


But wait, this isn't good enough for a few reasons:

- the rows have `ID` and `Parent` fields that are important!  But they're wrapped up in that `attributes` column.
- some rows have missing data that is represented by `.`.  But it should be parsed as a missing value.

What we need is a function:

```
def readGFF3File( filename ):
    """Read the GFF3 file with the given filename.  Return a dataframe
    wih ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns"
    result = # something
    # (populate result here)
    return result
```

This function higlights a few of our principles:

- it does one thing (keeping it simple)
- it is named for what it does.  (Function names should be verbs.)

**Question.** Can you write the above function?


**Note**. From a coding point of view it's a bit annoying that this function reads in a filename. Why should this code
be fussed with opening files? Also what if I wanted to send it data from somewhere else? We're doing it this way here
because that's what `pandas.read_table()` is expecting, but we will have to jump through some hoops to test it below.

**Note.** I tend to write functions this way a lot: specifically the result of the function is always called `result`,
the first line typically declares `result` and the last line is `return result`. The function's job is then to create
`result`.  This keeps things as simple as possible within the function.

Here's my attempt:

````
def readGFF3File( filename ):
    """Read the GFF3 file with the given filename.  Return a dataframe
    wih ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns"
    result = pandas.read_table(
        filename,
        comment = '#',
        names = [ 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes' ]
    )

    def parseAttributes( string ):
        """Parse a GFF3 attributes column (semi-colon separated key=value pairs)
        Return a dict of key/value pairs"""
        result = {}
        parts = string.split( ";" )
        for part in parts:
            keyValue = part.split( "=" )
            result[keyValue[0]] = keyValue[1]
        return result
    
    # do something here
    return result








Following 'keep it simple' let's focus on bits we can do first. It is a fixed format with multiple rows and 9 columns
that have names. This is what *data frames* were created to solve. In python, data frames are implemented in the
`pandas` package, so we can read in the data like this:

```
import pandas

```



