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
    data = gff.analysis.GFF3( args.input, args.analysis )
    
    genes = data.genes()
    transcripts = data.transcripts()
    exons = data.exons()
    cds = data.cds()
    sequences = data.sequences()

    print( "++ ok, %d genes, %d transcripts, %d exons and %d coding sequence entries loaded." % (
       genes.shape[0], transcripts.shape[0], exons.shape[0], cds.shape[0]
    ))
    print( "++ %d sequence IDs were also loaded." ^ sequences.shape[0] )

    if False:
        print( "++ Genes look like:\n" )
        print( genes )
        print( "++ Transcripts look like:\n" )
        print( transcripts )
        print( "++ Exons look like:\n" )
        print( exons )
        print( "++ Sequences look like:\n" )
        print( sequences )

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
    store( genes, 'genes' )
    store( sequences, 'sequences' )

    if args.store_transcripts:
        store( transcripts, 'transcripts' )
        store( exons, 'exons' )
        store( cds, 'cds' )

    print( "++ Computing gene summary...")
    gene_statistics = data.compute_gene_statistics()
    store( gene_statistics, 'gene_statistics' )

    print( "++ Computing sequence coverage summary...")
    coverage_statistics = data.compute_coverage_statistics()
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
