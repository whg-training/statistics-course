# gff.py
# This file implements the function parse_gff3_to_dataframe()
# and a number of helper functions.
import pandas

def gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""
    result = _read_gff3_using_pandas( file )
    extract_attributes_to_columns( result, ['ID', 'Parent', 'Name' ] )
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

def extract_attributes_to_columns( data, attributes_to_extract = [ 'ID', 'Parent' ], offset = 0 ):
    # The original function called parse_attributes twice - and was slow.
    # The second version called it once to unpack them all, but used masses of memory.
    # This third version uses regular expression to parse the attributes, get the value
    # and remove the extracted ones from the string directly - without duplicating the memory.
    import re
    def getAttribute( entry, regexp ):
        m = re.search( regexp, entry )
        return None if m is None else m.group(1)
    def removeAttribute( entry, regexp ):
        return re.sub( regexp, "", entry )

    for i in range( 0, len(attributes_to_extract)):
        attribute = attributes_to_extract[i]
        regexp = re.compile( "%s=([^;]+);?" % attribute )
        data.insert( i + offset, attribute, data['attributes'].apply( lambda entry: getAttribute( entry, regexp ) ))
        # Delete the field from the current attributes
        # (I could not get .transform() work here for some reason)
        data.loc[ :, 'attributes'] = data['attributes'].apply( lambda entry: removeAttribute( entry, regexp ))

def parse_sequences_from_gff_metadata( file ):
    """GFF3 files from the Ensembl ftp site list sequences and their lengths in the file metadata.
    This function parses this information and returns it as a pandas dataframe.
    It's use may be specific to the Ensembl files."""
    import pandas
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
    return pandas.DataFrame( result )

