# gff.py
# This file implements three functions:
# parse_gff3_to_dataframe() (the pandas version explored in the tutorial)
# parse_gff3_to_list() (a basic python version which does something similar)
# parse_attributes() which is a helper function to parse the attributes column.
# parse_sequences_from_gff_metadata() (reads sequence lengths from the metadata in Ensembl .gff files)

# I've used a particular style here:
# - names_are_written_using_underscores (asOpposedToCamelCase, which I could have done)
# - I've used nested functions to keep related bits of code together.-
# - all functions create a `result` variable on their first line and their job is to construct and return it.
#
# I've also written this in the top-down style, where the main function goes first followed by
# the functions it uses etc.
# 
def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""

    result = _read_gff3_using_pandas( file ) # this is defined below

    # The original function called parse_attributes twice - and was slow.
    # I got rid of this and now unpack them once:
    attributes = result['attributes'].apply( parse_attributes )
    # `attributes` is a list of dicts representing the unpacked attributes strings.
    
    # Extract desired columns
    attributes_to_extract = [ 'ID', 'Parent', 'Name', 'biotype' ]
    for attribute in attributes_to_extract:
        # Add a new column with the given name to the dataframe. The code on the
        # right is a python 'list comprehension', which is a quick way to write
        # something applied to every element of a list.  (So it is like pandas
        # .apply() )
        result[attribute] = [ entry.get( attribute, None ) for entry in attributes ]
        # One problem when working with column data formats is how to name
        # variables.  if the variable holding the column is called 'attributes',
        # like it is here, then what should I call the variable holding one
        # entry of the column? Here I'm calling it 'entry' but it's not ideal.

    def format_attributes( entry, exclude ):
        return ';'.join( "%s=%s" % ( key, value ) for key, value in entry.items() if key not in exclude )
    
    result['attributes'] = [ format_attributes( entry, exclude = attributes_to_extract ) for entry in attributes ]
    # reorder columns, because I want ID and Parent first
    result = result[ ['ID', 'Parent', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'Name', 'biotype', 'attributes'] ]

    return result

# functions starting with underscores are private to the file
def _read_gff3_using_pandas( file ):
    """Helper function to read the given GFF3 file into a dataframe, without any postprocessing."""
    import pandas
    result = pandas.read_table(
        file,
        comment = '#',
        names = [ 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes' ],
        na_values = ".",
        dtype = {
            'seqid': str,
            'source': str,
            'type': str,
            'start': int,
            'end': int,
            'score': float,
            'strand': str,
            'phase': str,
            'attributes': str
        }
    )
    return result

def parse_attributes( string ):
    """Helper function to parse a GFF3 attributes column (i.e. a string with semi-colon separated key=value pairs)
    and return a dict of the key/value pairs"""
    result = {}
    parts = string.split( ";" )
    for part in parts:
        keyValue = part.split( "=" )
        result[keyValue[0]] = keyValue[1]
    return result

def parse_sequences_from_gff_metadata( file ):
    """GFF3 files from the Ensembl ftp site list sequences and their lengths in the file metadata.
    This function parses this information and returns it as a pandas dataframe.
    It's use may be specific to the Ensembl files."""
    result = []
    for line in file:
        if line.startswith( '##sequence-region' ):
            parts = line.strip().split( " " )
            nameStartEnd = parts[-3:] # last 3 elements
            result.append({
                "seqid": nameStartEnd[0],
                "start": int( nameStartEnd[1] ),
                "end": int( nameStartEnd[2] )
            })
        elif not line[0] == '#':
            # quit when we meet the first non-metadata line
            break
    import pandas
    return pandas.DataFrame( result )
