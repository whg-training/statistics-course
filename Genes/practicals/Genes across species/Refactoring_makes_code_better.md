[Up to table of contents](README.md)

[Back to the previous page](Converting_gff_to_sqlite.md)

## Refactoring makes code better

How much do you like the code you wrote in `gff.py`?

Personally, I don't like [my original version](solutions/part1/gff_first_version.py) very much! Shorn of comments it
looks like this:

```
def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""

    def getID( attributes ):
        return parse_attributes( attributes ).get( 'ID', None )
    def getParent( attributes ):
        return parse_attributes( attributes ).get( 'Parent', None )
    result = _read_gff3_using_pandas( file )
    result['ID'] = result['attributes'].apply( getID )
    result['Parent'] = result['attributes'].apply( getParent )
    result = result[ ['ID', 'Parent', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'] ]
    return result
```

It's just... *complicated*.

The function should be simple.  It only does three things:

1. Loads the data
2. Extract the `ID` and `Parent` fields from the `attributes` column
3. Returns it.

So that is how we should write it:

```
def parse_gff3_to_dataframe( file ):
    result = read_gff3_using_pandas( file )
    add_ID_and_Parent( result )
    return result
```

That's more like it! This new version of the function is much shorter and simpler, and it's easy to understand. A good
sign is that (unlike the first version) it doesn't really need any comments because it's more or less obvious what it
does just from the function names.

Of course we have to write the `add_ID_and_Parent()` function, but [that's not hard](solutions/part1/gff.py). In the
process of writing this I realised that the [pandas `.insert()`
method](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.insert.html) can be used to put the columns in at
the start - so the last line of the original version (that reorders columns) isn't needed any more.

The `add_ID_and_Parent()` function follows a different pattern from the other functions we have written: it *mutates
its first argument*.  

Written this way, both the main function and the bits used for implementation are easy to test. Here is a [test using
python's unit test module]( (solutions/part1/test_gff.py) to show you what I mean.

**Note.** If you run the tests, you'll see that in fact my `parse_attributes()` function has a potential bug if it
encounters a badly-formatted attributes value that contains an equals sign in the value itself. Does this type of input
appear in the wild? Probably. (The [`partition()`](https://docs.python.org/3/library/stdtypes.html) function could be
used to fix this.)

**Challenge.** Refactor your code. Can you make it easy to understand?

Now [go back to the section on writing code](Getting_started_writing_some_code.md).
