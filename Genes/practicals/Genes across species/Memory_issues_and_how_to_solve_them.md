[Up to table of contents](README.md)
[Back to the previous page](Counting_genes_2.md)

## Memory issues - and how to solve them

So far in this tutorial we haven't thought much about memory (computer memory that is), or how best
we should use it. We've used python and pandas to implement our code, because that was easy and
gave us high-level features. And we haven't done anything special to keep memory usage low or work
on performance.

Unfortunately the data volumes used in many bioinformatics analyses are so large, this approach
can fall over quite quickly. And indeed if I use [version 2 of my code](solutions/part2/gff.py)] to load the gene
annotations:
```
import gff
data = gff.parse_gff3_to_dataframe("Homo_sapiens.GRCh38.104.chr.gff3" )
```

I find that the process is taking over 3.5Gb of RAM. (My laptop has 16Gb, so this still runs, but
if the data was any larger, it would start to fail.)

Moreover, seemingly innocuous changes to the code can make large differences to the memory usage.
My two versions of `gff.py` - [version 1](solutions/part1/gff.py) and [version
2](solutions/part2/gff.py) - do almost the same thing. The second one was tweaked to add a couple of
extra columns. But it turns out to use about 4 times the memory.

**Question.** How much memory does your version of the code use? (To find out, run the above in one
terminal, then in a second terminal run `top -u <username> -o RES` (linux) or `top -U <username> -o
mem` (Mac OS). You should see your program climbing the rankings as it loads the data.

### Some solutions

The main ways to solve this type of problem are:

1. work with data subsets (specific genome regions, specific record types, smaller organisms, or work in chunks).
2. use memory-efficient data structures.
3. use things like sqlite that allow efficient access to on-disk data.
4. write code more carefully to control memory usage.

In our problem we can exploit 1 and 3 by loading data directly from the sqlite database. (This of
course assumes your `gff_to_sqlite.py` program didn't run out of memory.)  That is we can write:

```
import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
data = pandas.read_sql( "SELECT * FROM gff_data WHERE type == 'gene'", db )
```

Because this loads only the 'genes' records, it uses much less memory.

Another thing we can do is try to improve the code. If you [look at the input
file](What_gene_annotation_data_looks_like.md) for a moment you'll realise that a huge proportion
of the actual data is in that `attributes` field. That's the one we have to work hardest with to
parse. Not surprisingly, what we do with this field really matters in terms of memory usage. The
key difference between my two versions of the `parse_gff3_to_dataframe()` function is that [the
second version](solutions/part2/gff.py) creates a variable that stores all the unpacked attributes:

```
    attributes = result['attributes'].apply( parse_attributes )
```
This one line turns out to use up about 2.5Gb of memory!

In [solutions/part2/gff_lowmem.py](solutions/part2/gff_lowmem.py) is one possible solution to this.
Instead of fully parsing the whole attribute string (which has a load of stuff we don't care
about), this version uses [regular expressions](https://en.wikipedia.org/wiki/Regular_expression)
(implemented in the [python re module](https://docs.python.org/3/library/re.html)) to parse out the
fields of interest. And then uses the same regular expression to remove the field from the
attributes. The crucial point here, however, is not the regular expressions themselves (which are
just a tool to do the parsing) but the fact that the data is never copied around. This version is
back down to < 2Gb of memory.

### Memory vs. performance

Funnily enough this version is also almost twice as fast as the version not using regular
expression! This is not really surprising either, because the regular expression library is highly
optimised; and our original `parse_attributes()` version was pretty simple and spent quite a bit of
time splitting up fields (and allocating memory for them) that were never used.

### Implications

Considerations like this mean we need to add a fourth item to our list of coding aims for
bioinformatics:

- it ought to work
- it ought to not take too long to do it
- it ought to be obvious what it does
- it ought to be frugal with memory

### Memory use challenge

In fact our program `gff_to_sqlite.py` does not need to use much memory at all. If you would like
an additional challenge, try this one:

**Advanced challenge:** write a version of `gff_to_sqlite.py` that has low memory usage.

**Hints**:

- The progam could process rows one at a time, instead of loading them all in at the top.

- The [the pure python version](gff_to_sqlite_python_version.py) is a good place to look for the
  needed sql statements.

- To make this slick, you should aim to commit the rows to the database in chunks (say of 10,000 rows). Then call
  `db.commit()` after each chunk to write the data.
  
- Don't forget to `db.commit()` after the last chunk
  
Is this version faster or slower than the original?
