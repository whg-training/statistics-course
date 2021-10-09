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
