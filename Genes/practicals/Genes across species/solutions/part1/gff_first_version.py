# gff.py
# This file implements the function parse_gff3_to_dataframe()
#
# and a number of helper functions.

# I've used a particular style here:
# - names_are_written_using_underscores (asOpposedToCamelCase)
# - I've generally made functions create a `result` variable, and their job is to build & return it.
# - I've also written it in the top-down style, in which the highest-level function goes first followed by the functions it uses etc.
# 
def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns"""

    # These are two helper functions to extract ID and Parent fields:
    def getID( attributes ):
        return parse_attributes( attributes ).get( 'ID', None )
    def getParent( attributes ):
        return parse_attributes( attributes ).get( 'Parent', None )

    result = read_gff3_using_pandas( file ) # this is defined below

    # Extract ID and Parent columns using the `apply()` dataframe method.
    result['ID'] = result['attributes'].apply( getID )
    result['Parent'] = result['attributes'].apply( getParent )

    # reorder columns, because I want ID and Parent first
    result = result[ ['ID', 'Parent', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'] ]

    return result

# functions starting with underscores are private to the file
def read_gff3_using_pandas( file ):
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

def parse_gff3_to_list( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a list of python `dicts` with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns."""
    result = []
    for line in file:
        if line[0] != '#':
            result.append( _from_gff3_line_to_dict( line ) )
    return result

def from_gff3_line_to_dict( line ):
    """Helper function to parse a single line of a GFF file into a python dict"""
    fields = line.strip().split( "\t" )
    assert len( fields ) == 9 # sanity check
    result = {
        "seqid": fields[0],
        "source": fields[1],
        "type": fields[2],
        "start": None if fields[3] == "." else int( fields[3] ),
        "end": None if fields[4] == "." else int( fields[4] ),
        "score": None if fields[5] == "." else float(fields[5]),
        "strand": None if fields[6] == "." else fields[6],
        "phase": None if fields[7] == "." else fields[7],
        "attributes": parse_attributes( fields[8] ),
        "attributes_string": fields[8] # keep this around as useful.
    }
    result["ID"] = result["attributes"].get( "ID", None )
    result["Parent"] = result["attributes"].get( "Parent", None )
    return result
