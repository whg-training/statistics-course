"""
gff_to_sqlite.py
This program reads data in GFF3 format and stores it in a table in a sqlite file, indexed by the ID.
Use it like this:

    python gff_to_sqlite.py --input <path to gff file> --output <path to output file> --analysis <name>

"""
import argparse, gff, sqlite3

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Convert a GFF3 file to sqlite3 format.
        The result will be a table with the GFF3 fields, and with ID, Parent, Name and buitype fields in columns.
        The resulting table will be indexed by the ID field for easy lookup."""
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
    return parser.parse_args()

def process( args ):
    print( "++ Loading gene annotation data from %s...\n" % args.input )
    gff3 = gff.parse.GFF3( args.input, args.analysis )
    
    genes = gff3.genes()
    transcripts = gff3.transcripts()
    exons = gff3.exons()
    cds = gff3.CDS()
    sequences = gff3.genome_sequences()

    print( "++ ok, %d genes, %d transcripts, %d exons and %d coding sequences loaded.\n" % (
       genes.shape[0], transcripts.shape[0], exons.shape[0], cds.shape[0]
    ))
    print( "++ Genes look like:\n" )
    print( genes )
    print( "++ Transcripts look like:\n" )
    print( transcripts )
    print( "++ Exons look like:\n" )
    print( exons )
    print( "++ Sequences look like:\n" )
    print( sequences )

    print( "++ Writing records to %s...\n" % args.output )
    db = sqlite3.connect( args.output )

    # In this version I have hard-coded `gff_data` and `sequences` table names.
    print( "++ Writing genes...\n" )
    genes.to_sql( 'genes', db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    print( "++ Writing transcripts...\n" )
    transcripts.to_sql( 'transcripts', db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    print( "++ Writing exons...\n" )
    exons.to_sql( 'exons', db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    print( "++ Writing cds...\n" )
    cds.to_sql( 'cds', db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    print( "++ Writing genome sequences...\n" )
    sequences.to_sql( 'sequences', db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    print( "++ Success.\n" )

# It is always good to let the user know what we are doing
# Both to confirm we're doing what they want, and because some of the steps
# can take some time (so we want to give them confidence it's still working).
# So we print out lots of messages!.print( "++ gff_to_sqlite_py" )
args = parse_arguments()
print( "++ Converting %s to %s...\n" % ( args.input, args.output ))
process( args )
print( "++ Thank you for using gff_to_sqlite.py" )
