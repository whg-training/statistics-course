"""
gff_to_sqlite.py
This program reads data in GFF3 format and stores it in a table in a sqlite file, indexed by the ID.
Use it like this:

    python gff_to_sqlite.py --input <path to gff file> --output <path to output file> --analysis <name>

"""
import argparse, gff_lowmem, sqlite3

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Convert a GFF3 file to sqlite3 format.
        The result will be a table with the GFF3 fields, and with ID and Parent fields in columns.
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
    print( "++ Loading genes data from %s...\n" % args.input )
    data = gff_lowmem.parse_gff3_to_dataframe( open( args.input ))
    print( "++ ok, %d records loaded, they look like:\n" % data.shape[0] )
    print( data )

    print( "++ Loading sequence lengths from %s...\n" % args.input )
    sequences = gff_lowmem.parse_sequences_from_gff_metadata( open( args.input ))
    print( "++ ok, %d records loaded, they look like:\n" % sequences.shape[0] )
    print( sequences )

    print( "++ Writing records to %s...\n" % args.output )
    db = sqlite3.connect( args.output )

    # In this version I have hard-coded `gff_data` and `sequences` table names.
    data.insert( 0, 'analysis', args.analysis )
    data.to_sql( "gff_data", db, index = False, if_exists = 'replace' if args.overwrite else 'append' )
    sequences.insert( 0, 'analysis', args.analysis )
    sequences.to_sql( "sequences", db, index = False, if_exists = 'replace' if args.overwrite else 'append' )

    print( "++ Success.\n" )

# It is always good to let the user know what we are doing
# Both to confirm we're doing what they want, and because some of the steps
# can take some time (so we want to give them confidence it's still working).
# So we print out lots of messages!.print( "++ gff_to_sqlite_py" )
args = parse_arguments()
print( "++ Converting %s to %s:%s...\n" % ( args.input, args.output, args.table ))
process( args )
print( "++ Thank you for using gff_to_sqlite.py" )
