## Memory issues and how to solve them
What's more, with these data volumes, seemingly innocuous changes to this type of code can start to make a huge
difference. My two versions of `gff.py` - [version 1](solutions/part1/gff.py) and [version 2](solutions/part2/gff.py) -
do almost the same thing. The second one just adds a couple of columns. But it turns out to use about 4 times the
memory.

Not surprisingly, in our program this is all to do with attribute parsing. They key difference between the two versions
is that the original version of `parse_gff3_to_dataframe()` **never stored the parsed/unpacked attributes strings**.
(Instead - [as profiling revealed](Converting_gff_to_sqlite.md) - it was wasting time by parsing them twice). But the
second version does do this: on line 25 it says:

```
    attributes = result['attributes'].apply( parse_attributes )
```
This one line turns out to use up about 2.5Gb of memory!

Considerations like this mean we need to add an additional aim when we are coding:

- it ought to work
- it ought to not take too long to do it
- it ought to be obvious what it does
- **it ought not to waste memory**

**Question.** How much memory is your version of the code using?

The main approaches to solve this problem are:

- work with data subsets (genome regions, specific record types, smaller organisms)
- use memory-efficient data structures (like the sqlite database we are using)
- write code more carefully to control memory usage

**Advanced challenge:** write a version of `gff_to_sqlite.py` that has low memory usage.

**Hints**:

- There's really no need for `gff_to_sqlite.py` to use any memory at all! It could process rows one at a time instead
  of loading them all in.

- The [the pure python version](gff_to_sqlite_python_version.py) is a good place to look for the needed sql statements.

- To make this slick, you should aim to commit the rows to the database in chunks (say of 10,000 rows). Then call
  `db.commit()` to write the data to the file after every chunk.

- Don't forget to commit the last chunk.
