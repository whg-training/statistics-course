def summarise_transcripts_and_exons_per_gene( genes, transcripts, exons ):
    """Given data frames of genes, transcripts and exons, return a dataframe with one row per gene
    and ID, number_of_transcripts, and average_number_of_exons columns."""
    transcript_summary = count_exons_per_transcript( transcripts, exons )
    return summarise_transcripts_per_gene( genes, transcript_summary )


def count_exons_per_transcript( transcripts, exons ):
    """Given data frames of transcripts and exons, return a dataframe with one row per transcripts
    and ID, Parent, and number_of_exons columns."""
    import pandas
    summary = _join_dataframes_by_ID_and_Parent( transcripts, exons, [ "ID", "Parent" ] )
    summary.rename(
        columns = { "ID_x": "ID", "Parent_x": "Parent", "ID_y": "exon_ID", "Parent_y": "exon_Parent" },
        inplace = True
    )
    result = summary.groupby( [ 'ID', 'Parent'], as_index = False ).agg(
        number_of_exons = pandas.NamedAgg(
            column = "exon_Parent",
            aggfunc = lambda x: x.notnull().sum() )
    )
    return result

def summarise_transcripts_per_gene( genes, transcript_summary ):
    """Given data frames of genes and a transcript-exon summary as returned by
    count_exons_per_transcript(), return a dataframe with one row per gene
    and ID, number_of_transcripts, and average_number_of_exons columns."""
    import pandas, numpy
    summary = _join_dataframes_by_ID_and_Parent( genes, transcript_summary, ["ID", "Parent", "number_of_exons"] )
    summary.rename( columns = { "ID_x": "ID", "Parent_x": "Parent", "ID_y": "transcript_ID", "Parent_y": "transcript_Parent" }, inplace = True )
    result = summary.groupby( ["analysis", "ID", "biotype", "Name", "seqid", "start", "end", "strand" ], dropna = False ).agg(
        number_of_transcripts = pandas.NamedAgg(
            column = "transcript_Parent",
            aggfunc = lambda x: x.notnull().sum()
        ),
        average_number_of_exons = pandas.NamedAgg(
            column = "number_of_exons",
            aggfunc = lambda x: numpy.mean(x)
        )
    )
    result.reset_index( inplace = True )
    return result

def _join_dataframes_by_ID_and_Parent( left, right, right_columns ):
    import pandas
    return pandas.merge(
        left,
        right[ right_columns ], 
           how = "outer",
           left_on = "ID",
           right_on = "Parent"
        )

def count_single_exon_genes( gene_summary ):
    import pandas
    result = gene_summary[ gene_summary['biotype'] == 'protein_coding' ].groupby( 'analysis' ).agg(
        total_genes = pandas.NamedAgg(
            column = 'ID',
            aggfunc = lambda x: x.notnull().sum()
        ),
        total_single_exon = pandas.NamedAgg(
            column = 'average_number_of_exons',
            aggfunc = lambda x: ( x == 1 ).sum()
        ),
        proportion = pandas.NamedAgg(
            column = 'average_number_of_exons',
            aggfunc = lambda x: (x == 1).sum() / x.notnull().sum()
        )
    )
    return result

