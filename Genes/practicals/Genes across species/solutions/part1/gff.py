# gff.py
# This file implements the function parse_gff3_to_dataframe()
# and a number of helper functions.
# This is the 'refactored' version in which the main function is only 3 lines long.

def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""
    result = read_gff3_using_pandas( file )
    add_ID_and_Parent( result )
    return result

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

def add_ID_and_Parent( data ):
    # These are two helper functions to extract ID and Parent fields:
    def getID( attributes ):
        return parse_attributes( attributes ).get( 'ID', None )
    def getParent( attributes ):
        return parse_attributes( attributes ).get( 'Parent', None )
    data.insert( 0, 'ID', data['attributes'].apply( getID ) )
    data.insert( 1, 'Parent', data['attributes'].apply( getParent ) )

def parse_attributes( string ):
    """Helper function to parse a GFF3 attributes column (i.e. a string with semi-colon separated key=value pairs)
    and return a dict of the key/value pairs"""
    result = {}
    parts = string.split( ";" )
    for part in parts:
        keyValue = part.split( "=" )
        result[keyValue[0]] = keyValue[1]
    return result

