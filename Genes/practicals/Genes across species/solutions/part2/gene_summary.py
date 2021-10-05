#
# This code is a python version of the SQL code to join
# exons to transcripts and transcripts to exons that we saw in the tutorial.
#
# The syntax  form the joins using the slightly complicated
# pandas syntax for groupby() operations with named columns.
# It is described here: https://pandas.pydata.org/docs/user_guide/groupby.html.
#

import pandas, sqlite3
db = sqlite3.connect( "genes.sqlite" )
genes = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'gene' )", db )
transcripts = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'mRNA' )", db )
exons = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand FROM gff_data WHERE type IN ( 'exon' )", db )

def count_exons_per_transcript( transcripts, exons ):
    summary = join_dataframes_by_ID_and_Parent( transcripts, exons, [ "ID", "Parent" ] )
    summary.rename( columns = { "ID_x": "ID", "Parent_x": "Parent", "ID_y": "exon_ID", "Parent_y": "exon_Parent" }, inplace = True )
    result = summary.groupby( [ 'analysis', 'ID', 'Parent'], as_index = False ).agg(
        number_of_exons = pandas.NamedAgg( column = "exon_ID", aggfunc = numpy.size )
    )
    return result

def join_dataframes_by_ID_and_Parent( left, right, right_columns ):
    return pandas.merge(
        left,
        right[ right_columns ], 
           how = "outer",
           left_on = "ID",
           right_on = "Parent"
        )

def summarise_transcripts_per_gene( genes, transcript_summary ):
    summary = join_dataframes_by_ID_and_Parent( genes, transcript_summary, ["ID", "Parent", "number_of_exons"] )
    summary.rename( columns = { "analysis_x": "analysis", "ID_x": "ID", "Parent_x": "Parent", "ID_y": "transcript_ID", "Parent_y": "transcript_Parent" }, inplace = True )
    
    result = summary.groupby( [ "analysis", "ID" ] ).agg(
        number_of_transcripts = pandas.NamedAgg( column = "ID", aggfunc = numpy.size ),
        average_number_of_exons = pandas.NamedAgg( column = "number_of_exons", aggfunc = numpy.mean ),
    )
    result.reset_index( inplace = True )
    return result

# Note: if we have got this right, the results should be identical to the SQL
# gene_summary_view developed in the tutorial.  Is it?

