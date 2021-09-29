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

    # These are two helper functions to extract ID and Parent fields:
    def getID( attributes ):
        return parse_attributes( attributes ).get( 'ID', None )
    def getParent( attributes ):
        return parse_attributes( attributes ).get( 'Parent', None )

    result = _read_gff3_using_pandas( file ) # this is defined below

    # Extract ID and Parent columns using the `apply()` dataframe method.
    result['ID'] = result['attributes'].apply( getID )
    result['Parent'] = result['attributes'].apply( getParent )

    # reorder columns, because I want ID and Parent first
    result = result[ ['ID', 'Parent', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'] ]

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

def parse_gff3_to_list( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a list of python `dicts` with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns."""
    result = []
    for line in file:
        if line[0] != '#':
            result.append( _from_gff3_line_to_dict( line ) )
    return result

def _from_gff3_line_to_dict( line ):
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
