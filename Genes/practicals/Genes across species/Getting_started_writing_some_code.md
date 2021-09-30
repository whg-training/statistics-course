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

**Use a package.** The [pandas library](https://pandas.pydata.org) is an obvious one to try here. The GFF3 data is
basically a tabular data format (many rows x 9 named columns, at least if we don't unpack `attributes`) and that's a
good fit for a dataframe (which is what pandas provides).  

**Use raw python.** This works fine as well, and it has some advantages in flexibility.

In this tutorial we'll focus on the pandas version, as it gives us lots of tools for analysis. We'll develop a little
python module that carries out the task.

## Diving straight in - parsing data

If you [looked at the gene annotation data](What_gene_annotation_data_looks_like.md), you'll know that it comes in rows
of data that are tab-delimited, but that it is also relational (meaning that the records refer to each other, via the
`Parent` attribute. Moreover since exons are associated with transcripts, which are in turn associated with genes,
we'll have to build some form of data hierarchical data structure for this.

That sounds complex, but we can break off a manageable bit of the job by just focussing on getting the data in (c.f.
"keep it simple"). So let's do the simplest thing possible and start by writing a function:

```
def parse_gff3_to_dataframe( data ):
    """Parse data in GFF3 format, and return a pandas dataframe"""
    result = (do some work here)
    # more code here to add to result
    return (result )
```

**Note:** If you are not familiar with python syntax, now would be a good time to refresh via any of the available
tutorials. The above code defines a function, and shows a documentation comment (the `"""..."""` bit) and also has code
comments (starting with `#`).

The function above already illustrates a couple of things that might be helpful if you're not used to writing code.
First, the name is very clear what this function does - indeed the documentation comment is pretty useless at the
moment.  Second, even before we've written it, the function follows a very simple pattern: it creates a new thing (the
result of the function, so it is called `result`) and the last line returns it. All the function has to do is build
`result` - simple!

The other thing is that this function is already reasonably testable.  Look, here is a test:

```
test_data = """##gff-version 3
#description: test data
chr1\tme\tgene\t1\t1000\t.\t+\t.\tID=gene1;other_data=stuff
chr1\tme\texon\t10\t900\t.\t+\t.\tID=gene1.1;Parent=gene1
"""

import io, math

# run the function:
data = parse_gff3_to_dataframe( io.StringIO( test_data ))

# test it:
assert data['seqid'][0] == 'chr1'
assert data['strand'][0] == '+'
assert data['attributes'][0] == 'ID=gene1;other_data=stuff'

assert data['start'][1] == 10 # an integer

assert math.isnan( data['score'][1] ) # Because a "." should be missing data in GFF spec

assert data['ID'][0] == 'gene1'
assert data['ID'][1] == 'gene1.1'
assert data['Parent'][1] == 'gene1'
# etc.
```

(**Note:** In an ideal world we could pass in the data directly. However it is a bit annoying to write a function that
works both with a file and a string. The `io.StringIO()` bit above just turns the input data into a file-like object -
you can pretty much ignore it.

This leads us to a:

**Challenge:** write `parse_gff3_to_dataframe()`. Can you write this so all the tests pass?  Or at least some of them?

**Hints:** This test does have a few complexities to it.  Here are some points to think about that might be helpful:

- pandas has a [`read_table` function](https://pandas.pydata.org/docs/reference/api/pandas.read_table.html) that is a
  good way to get the data in.  It has many arguments that control how the data gets in there.

- The data itself doesn't have column names in. But the test requires them - you have to get them in somehow.

- To pass the tests you have to pay attention to missing data! Check the [GFF
  spec](https://m.ensembl.org/info/website/upload/gff3.html) for how these are represnted.

- Also, some of the columns [have different data types](https://m.ensembl.org/info/website/upload/gff3.html) - to pass
  the tests you also have to get the types right.

- The last three rows of the test use the `ID` and `Parent` fields. If you [examine the data]() you'll see these aren't
  columns in GFF3 but live inside the complicated `attributes` column. The test above is asking that you pull these out
  into new columns.
  
**More hints:**

- Be careful about whitespace. Python really cares about your indentation. In particular, you should decide at the
  start if you are going to indent with spaces or tabs and stick with it (otherwise you will get all sorts of strange syntax
  errors).

- It is very worth looking around for an editor you like. On the Mac, [Textmate](https://macromates.com) or
  [Sublime](https://www.sublimetext.com) or [github Atom](https://atom.io) might be good choices. On Windows I'm less
  sure - [Notepad++](https://notepad-plus-plus.org) or [Atom](https://atom.io) might work. Sublime and Atom also work
  on linux.  Of course there are also many others you can use.

- It's a good idea to write little helper functions to do bits of the task. (For example, it might be useful to write
 a `parse_attributes()` function that parses a semi-colon-separated list of key=value pairs and return a dict, e.g. like this:
```
> parse_attributes( 'ID=gene1;other_data=stuff' )
{
    "ID": "gene1",
    "other_data": "stuff"
}
```

- `.split()` can be used to split strings.  For example `"Hello;world".split( ";" ) == [ "hello", "world" ]``

- If you want to apply a function to every row of a data frame - the [pandas `.apply()`
  method](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.apply.html) is the ticket.
  For example, `int` is a function, so you could do:
```
    df['start'].apply( int )
```

- How do handle missing data?  One way is to use the python [syntax for conditional expressions](https://mail.python.org/pipermail/python-dev/2005-September/056846.html):
```
# convert value to an int
value = None if value == "." else int(value)
```
(But for this task an easier way may be to exploit the `na_values` argument of `pandas.read_table`.)

**Yet more hints:** The `solutions/part1/gff.py` file contains my solution to this. Feel free to have a look / steal code. As
a comparison it also implements a similar pure python version called `parse_gff3_to_list()`. (There are lots of other
ways to write this - for ecample, gathering the code into a class might be sensible - but I've gone with functions for
simplicity.)

### A note on writing tests first

I'd hazard a guess that not many people writing scientific code (including me) actually write their tests first in the
way we did above. It is however a very useful approach, mainly because it forces you to think about how your code will
be used before you spend the effort of writing it. So by definition it tends to produce functions etc. that are easy to
use.

Also, if you do this you'll find all your code is tested - for free! (If this sounds like I'm convincing myself, it's
because I hardly ever do this - but I wish I did.)

### Using the code

Do you have a working function `parse_gff3_to_dataframe()`?  Then let's [turn this into a useful program](Converting_gff_to_sqlite.md).
