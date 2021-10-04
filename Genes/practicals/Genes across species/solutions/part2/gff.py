# gff.py
# This file implements the function parse_gff3_to_dataframe()
# and a number of helper functions.

def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""
    result = read_gff3_using_pandas( file )
    extract_attributes_to_columns( result, ['ID', 'Parent', 'Name', 'biotype'] )
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

def extract_attributes_to_columns( data, attributes_to_extract = [ 'ID', 'Parent' ] ):
    # The original `add_ID_and_Parent()` function ultimately called parse_attributes twice - and was slow.
    # To fix this I now unpack them once at the top:
    attributes = data['attributes'].apply( parse_attributes )
    # `attributes` is now a list of dicts representing the unpacked attributes strings.

    for i in range( 0, len( attributes_to_extract )):
        attribute = attributes_to_extract[i]
        # The code in the 3rd argument below is a python 'comprehension'.  This is
        # a quick way to write something applied to every element of a list - in this case, to
        # extract the relevant attribute from the unpacked attributes.
        data.insert( i, attribute, entry.get( attribute, None ) for entry in attributes )

    # Now replace the original column with the remaining attributes:
    def format_attributes( entry, exclude ):
        return ';'.join( "%s=%s" % ( key, value ) for key, value in entry.items() if key not in exclude )
    result['attributes'] = [ format_attributes( entry, exclude = attributes_to_extract ) for entry in attributes ]

def parse_attributes( string ):
    """Helper function to parse a GFF3 attributes column (i.e. a string with semi-colon separated key=value pairs)
    and return a dict of the key/value pairs"""
    result = {}
    parts = string.split( ";" )
    for part in parts:
        keyValue = part.split( "=" )
        result[keyValue[0]] = keyValue[1]
    return result
