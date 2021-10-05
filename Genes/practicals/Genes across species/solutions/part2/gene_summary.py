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

transcript_summary = pandas.merge(
   transcripts,
   exons[["ID", "Parent"]], 
   how = "outer",
   left_on = "ID",
   right_on = "Parent"
).groupby( ['ID_x', 'Parent_x'], as_index = False ).agg(
    number_of_exons = pandas.NamedAgg( column = "ID_y", aggfunc = numpy.size )
)

gene_summary = pandas.merge(
    genes,
    transcript_summary,
    how = "outer",
    left_on = "ID",
    right_on = "Parent_x"
).groupby( [ "analysis", "ID" ] ).agg(
    number_of_transcripts = pandas.NamedAgg( column = "ID", aggfunc = numpy.size ),
    average_number_of_exons = pandas.NamedAgg( column = "number_of_exons", aggfunc = numpy.mean ),
)

# The table above is identical to the SQL gene_summary_view developed in the tutorial.
