# gff.py
# This file implements the function parse_gff3_to_dataframe()
# and a number of helper functions.
import pandas

def parse_gff3_to_dataframe( file ):
    """Read GFF3-formatted data in the specified file (or file-like object)
    Return a pandas dataframe with ID, Parent, seqid, source, type, start, end, score, strand, phase, and attributes columns.
    The ID and Parent are extracted from the attributes columns, and the dataframe is indexed by ID"""
    result = read_gff3_using_pandas( file )
    extract_attributes_to_columns( result, ['ID', 'Parent', 'Name', 'biotype'] )
    return result

# functions starting with underscores are private to the file
def read_gff3_using_pandas( file ):
    """Helper function to read the given GFF3 file into a dataframe, without any postprocessing."""
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
        data.insert( i, attribute, [entry.get( attribute, None ) for entry in attributes ])

    # Now replace the original column with the remaining attributes:
    def format_attributes( entry, exclude ):
        return ';'.join( "%s=%s" % ( key, value ) for key, value in entry.items() if key not in exclude )
    data['attributes'] = [ format_attributes( entry, exclude = attributes_to_extract ) for entry in attributes ]

def parse_attributes( string ):
    """Helper function to parse a GFF3 attributes column (i.e. a string with semi-colon separated key=value pairs)
    and return a dict of the key/value pairs"""
    result = {}
    parts = string.split( ";" )
    for part in parts:
        keyValue = part.split( "=" )
        result[keyValue[0]] = keyValue[1]
    return result

def summarise_genes( genes, transcripts, exons ):
    """Given data frames of genes, transcripts and exons, return a dataframe with one row per gene
    and ID, number_of_transcripts, and average_number_of_exons columns."""
    transcript_summary = count_exons_per_transcript( transcripts, exons )
    return summarise_transcripts_per_gene( genes, transcript_summary )


def count_exons_per_transcript( transcripts, exons ):
    """Given data frames of transcripts and exons, return a dataframe with one row per transcripts
    and ID, Parent, and number_of_exons columns."""
    summary = join_dataframes_by_ID_and_Parent( transcripts, exons, [ "ID", "Parent" ] )
    summary.rename(
        columns = { "ID_x": "ID", "Parent_x": "Parent", "ID_y": "exon_ID", "Parent_y": "exon_Parent" },
        inplace = True
    )
    result = summary.groupby( [ 'ID', 'Parent'], as_index = False ).agg(
        number_of_exons = pandas.NamedAgg(
            column = "exon_ID",
            aggfunc = lambda x: x.notnull().sum() )
    )
    return result

def summarise_transcripts_per_gene( genes, transcript_summary ):
    """Given data frames of genes and a transcript-exon summary as returned by
    count_exons_per_transcript(), return a dataframe with one row per gene
    and ID, number_of_transcripts, and average_number_of_exons columns."""
    import numpy
    summary = join_dataframes_by_ID_and_Parent( genes, transcript_summary, ["ID", "Parent", "number_of_exons"] )
    summary.rename( columns = { "ID_x": "ID", "Parent_x": "Parent", "ID_y": "transcript_ID", "Parent_y": "transcript_Parent" }, inplace = True )
    result = summary.groupby( [ "ID" ] ).agg(
        number_of_transcripts = pandas.NamedAgg(
            column = "transcript_ID",
            aggfunc = lambda x: x.notnull().sum()
        ),
        average_number_of_exons = pandas.NamedAgg(
            column = "number_of_exons",
            aggfunc = lambda x: numpy.mean(x)
        )
    )
    result.reset_index( inplace = True )
    return result

def join_dataframes_by_ID_and_Parent( left, right, right_columns ):
    return pandas.merge(
        left,
        right[ right_columns ], 
           how = "outer",
           left_on = "ID",
           right_on = "Parent"
        )

def summarise_genes_python_version( genes, transcripts, exons ):
    # This is a python-based implementation of the gene summary code.
    # It works by building an array data structure that keeps
    # track of exons per transcript, and then transcripts per gene.
    # The results are returned as a pandas data frame.

    # Here I am using nested functions to keep everything together:
    def count_exons_per_transcript( transcripts, exons ):
        result = [ None ] * transcripts.shape[0]
        transcript_ids = {}
        for i, transcript in transcripts.iterrows():
            transcript_ids[ transcript['ID'] ] = i
            result[i] = {
                "ID": transcript['ID'],
                'Parent': transcript['Parent'],
                "exons": []
            }
        for i, exon in exons.iterrows():
            transcript_i = transcript_ids[ exon['Parent'] ]
            result[transcript_i]['exons'].append( exon['ID'] )
        return result

    def summarise_transcripts_per_gene( genes, transcript_summary ):
        result = [ None ] * genes.shape[0]
        gene_ids = {}
        for i, gene in genes.iterrows():
            gene_ids[ gene['ID'] ] = i
            result[i] = {
                "ID": gene['ID'],
                "transcripts": []
            }
        for transcript in transcript_summary:
            gene_id = gene_ids[ transcript['Parent'] ]
            result[ gene_id ]['transcripts'].append( transcript )
        return result

    def compute_average_number_of_exons( gene ):
        transcripts = gene['transcripts']
        if len( transcripts ) == 0:
            return None
        count = sum( len( transcript['exons'] ) for transcript in transcripts )
        return count / len( transcripts )

    transcript_summary = count_exons_per_transcript( transcripts, exons )
    gene_summary = summarise_transcripts_per_gene( genes, transcript_summary )
    result = [
        {
            "ID": gene['ID'],
            "number_of_transcripts": len( gene['transcripts'] ),
            "average_number_of_exons": compute_average_number_of_exons( gene )
        }
        for gene in gene_summary
    ]
    return pandas.DataFrame( result )

