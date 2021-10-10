#
# This code is a python version of the SQL code to join
# exons to transcripts and transcripts to exons that we saw in the tutorial.
#
# The syntax  form the joins using the slightly complicated
# pandas syntax for groupby() operations with named columns.
# It is described here: https://pandas.pydata.org/docs/user_guide/groupby.html.
#

import pandas, sqlite3, gff
db = sqlite3.connect( "genes.sqlite" )
genes = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'gene' )", db )
transcripts = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand, Name, biotype FROM gff_data WHERE type IN ( 'mRNA' )", db )
exons = pandas.read_sql( "SELECT analysis, ID, Parent, seqid, start, end, strand FROM gff_data WHERE type IN ( 'exon' )", db )

# This takes a while...
gene_summary = gff.summarise_genes( genes, transcripts, exons )

# Add gene length because we didn't do that earlier.
gene_summary['length'] = gene_summary['end'] - gene_summary['start']

# Find the extremes
extremes = gene_summary.groupby( 'analysis' ).agg(
    max_size = pandas.NamedAgg(
        column = "length",
        aggfunc = lambda x: x.max()
    ),
    min_size = pandas.NamedAgg(
        column = "length",
        aggfunc = lambda x: x.min()
    ),
    max_exons = pandas.NamedAgg(
        column = "average_number_of_exons",
        aggfunc = lambda x: x.max()
    ),
    max_transcripts = pandas.NamedAgg(
        column = "number_of_transcripts",
        aggfunc = lambda x: x.max()
    )
)

# Find genes with highest number of transcripts...
pandas.merge(
    left = extremes,
    right = gene_summary,
    left_on = [ 'analysis', 'max_transcripts' ],
    right_on = [ 'analysis', 'number_of_transcripts' ]
)

# ...or highest length...
pandas.merge(
    left = extremes,
    right = gene_summary,
    left_on = [ 'analysis', 'max_length' ],
    right_on = [ 'analysis', 'length' ]
)

# ...or highest number of exons...
pandas.merge(
    left = extremes,
    right = gene_summary,
    left_on = [ 'analysis', 'max_exons' ],
    right_on = [ 'analysis', 'average_number_of_exons' ]
)

# Find the proportion of single-exon genes
# We use the average exon count for simplicity - this should really be improved.
single_exon_count = gene_summary[ gene_summary['biotype'] == 'protein_coding' ].groupby( 'analysis' ).agg(
    single_exon_count = pandas.NamedAgg(
        column = 'average_number_of_exons',
        aggfunc = lambda x: ( x == 1 ).sum()
    ),
    total = pandas.NamedAgg(
        column = 'average_number_of_exons',
        aggfunc = lambda x: x.notnull().sum()
    ),
    proportion = pandas.NamedAgg(
        column = 'average_number_of_exons',
        aggfunc = lambda x: (x == 1).sum() / x.notnull().sum()
    )
)

# make plots here?
