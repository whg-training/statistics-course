def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""

    result = _read_gff3_using_pandas( file ) # this is defined below

    # The original function called parse_attributes twice - and was slow.
    # The second version called it once to unpack them all, but used masses of memory
    # This third version uses regular expression to parse the attributes, get the value
    # and remove the extracted ones from the string directly - without duplicating the memory.
    import re # for regular expressions
    def getAttribute( entry, regexp ):
        m = re.search( regexp, entry )
        return None if m is None else m.group(1)
    def removeAttribute( entry, regexp ):
        return re.sub( regexp, "", entry )

    print( "Processing attributes..." )
    for name in [ 'ID', 'Parent', 'Name', 'biotype' ]:
        regexp = re.compile( "%s=([^;]+);?" % name )
        result[name] = result['attributes'].apply( lambda entry: getAttribute( entry, regexp ) )
        # Delete the field from the current attributes
        # (I could not .transform() work here for some reason)
        result['attributes'] = result['attributes'].apply( lambda entry: removeAttribute( entry, regexp ))
    
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
