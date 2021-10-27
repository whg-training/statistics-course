"""
gff_to_sqlite.py
This program reads data in GFF3 format and gathers summary statistics for genes.  Results are stored in a sqlite file.
Use it like this:

    python summarise_gff3.py --input <path to gff file> --output <path to output file> --analysis <name>

"""
import argparse, gff, sqlite3

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Gather summary statistics for genes from a GFF3 file to sqlite3 format.
        The result, including a table of genes, will be stored in a sqlite file for downstream analysis."""
    )
    parser.add_argument(
        '--analysis',
        help ='A name to give this analysis, e.g. species name',
        required = True
    )
    parser.add_argument(
        '--input',
        help ='The path of a GFF3-formatted file to work with',
        required = True
    )
    parser.add_argument(
        '--output',
        help ='The path of a sqlite3 file to output results to.',
        required = True
    )
    parser.add_argument(
        '--overwrite',
        action = "store_true",
        help ='If specified, overwrite the table with this data.  Otherwise data will be appended.'
    )
    parser.add_argument(
        '--store_transcripts',
        action = "store_true",
        help = """If specified, store transcripts, exons and coding sequence records in the output database.
        This can be useful but takes lots of space in the file."""
    )
    return parser.parse_args()

def process( args ):
    print( "++ Loading gene annotation data from %s...\n" % args.input )
    data = gff.analysis.UnpackedGFF.from_gff3( args.input, args.analysis )
    print( data.summary() )

    print( "++ Writing records to %s..." % args.output )
    db = sqlite3.connect( args.output )

    print( "++ Writing gene summary..." )
    def store( dataframe, table_name ):
        print( "++ Writing %s..." % table_name, end = '' )
        dataframe.to_sql(
            table_name,
            db,
            index = False,
            if_exists = 'replace' if args.overwrite else 'append'
        )
        print( " Ok, %d records written." % dataframe.shape[0] )
    # We just write the summary which has all the relevant fields
    store( data.genes, 'genes' )
    store( data.sequences, 'sequences' )

    if args.store_transcripts:
        store( data.transcripts, 'transcripts' )
        store( data.exons, 'exons' )
        store( data.cds, 'cds' )

    print( "++ Summarising genes, transcripts, and exons...")
    gene_summary = gff.analysis.summarise_genes( data )
    print( gene_summary['statistics'] )
    store( gene_summary['per_gene'], 'gene_summary' )
    store( gene_summary['statistics'], 'gene_statistics' )

    print( "++ Summarising gene coverage...")
    coverage_statistics = gff.analysis.compute_coverage_statistics( data )
    store( coverage_statistics, 'coverage_statistics' )

    print( "++ Success.\n" )

# It is always good to let the user know what we are doing
# Both to confirm we're doing what they want, and because some of the steps
# can take some time (so we want to give them confidence it's still working).
# So we print out lots of messages!.print( "++ summarise_gff3_py" )
args = parse_arguments()
print( "++ Converting %s to %s...\n" % ( args.input, args.output ))
process( args )
print( "++ Thank you for using summarise_gff3.py" )
