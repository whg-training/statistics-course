## Writing a useful conversion program.

In the [first part of this tutorial](Getting_started_writing_some_code.md) you will have written some code to parse a
GFF3 file and output a pandas dataframe.  I'll assume your function is called `parse_gff3_to_dataframe()`.

Your function is already useful! To demonstrate this, let's convert the GFF file into a different format - a [sqlite
database](https://www.sqlite.org) database. `sqlite` is a very useful file format that works like a full database, but
doesn't require a server or anything like that - it's just a single file on the filesystem. But unlike a flat file, it
is relational, making it easy to link records together. (This is [similar to pandas
functionality]([https://pandas.pydata.org/docs/getting_started/comparison/comparison_with_sql.html#compare-with-sql-join
) but doesn't need to load all the data into memory.)

To get started, put your function in a file called `gff.py`.  If you start python you should now be able to write:
```
import gff
data = gff.parse_gff3_into_dataframe('PlasmoDB-54_Pfalciparum3D7.head.gff')
print(data)
```

(Congratulations - you've written a python module!)

### Writing a command-line program

We will use your  write a program called `gff_to_sqlite_dataframe_version.py` which

- uses the [argparse](https://docs.python.org/3/library/argparse.html) library to process command-line arguments for the input and output files.
- loads the data using your function above.
- stores the data in the output sqlite file.

This turns out to be pretty easy. Start by importing the relevant stuff, and writing a function to parse the arguments.
The [argparse documentation](https://docs.python.org/3/library/argparse.html) is a bit confusing I find, so here is a
quick head start:

```
import argparse, sqlite3
import gff

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Convert a GFF3 file to sqlite3 format.
        The result will be a table with the GFF3 fields, and with ID and Parent fields in columns.
        The resulting table will be indexed by the ID field for easy lookup."""
    )
    parser.add_argument(
        '--input',
        help ='The path of a GFF3-formatted file to work with',
        required = True
    )

    # add other needed arguments here

    return parser.parse_args()

args = parse_arguments()
```

If you put this in a file called `gff_to_sqlite.py` then you can run it like this:
```
````

