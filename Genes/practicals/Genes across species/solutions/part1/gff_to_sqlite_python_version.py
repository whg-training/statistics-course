"""
gff_to_sqlite.py
This program reads data in GFF3 format and stores it in a table in a sqlite file, indexed by the ID.
Use it like this:

    python gff_to_sqlite_python_version.py --input <path to gff file> --output <path to output file> --analysis <name>

"""
import argparse, gff, sqlite3

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Convert a GFF3 file to sqlite3 format.
        The result will be a table with the GFF3 fields, and with ID and Parent fields in columns.
        The resulting table will be indexed by the ID field for easy lookup."""
    )
    parser.add_argument(
        '--analysis',
        help ='A name to give this analysis.  (I find this useful to avoid losing track of results).',
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
    data = gff.parse_gff3_to_list( open( args.input ))
    print( "++ ok, %d records loaded, first few look like:\n" % len(data))
    print( data[0:min(10,len(data))] )

    print( "++ Writing records to %s:%s...\n" % ( args.output, args.table ))
    db = sqlite3.connect( args.output )

    if args.overwrite:
        db.execute( "DROP TABLE IF EXISTS `%s`" % args.table )

    db.execute( """CREATE TABLE IF NOT EXISTS `%s` (
        `analysis` TEXT NOT NULL,
        `ID` TEXT,
        `Parent` TEXT,
        `seqid` TEXT,
        `source` TEXT,
        `type` TEXT,
        `start` INTEGER,
        `end` INTEGER,
        `score` REAL,
        `strand` TEXT,
        `phase` TEXT,
        `attributes` TEXT
    ) ;""" % args.table )
    
    db.executemany(
        "INSERT INTO `%s` VALUES( '%s', :ID, :Parent, :seqid, :source, :type, :start, :end, :score, :strand, :phase, :attributes_string )" % ( args.table, args.analysis ),
        data
    )
    db.commit()
    print( "++ Indexing...\n" )
    #print( "++ Indexing ID field...\n" )
    # Indexing is a good idea.  But if we are doing multiple datasets
    # it is faster to index them after, by running this SQL on the db:
    #db.execute( "CREATE INDEX IF NOT EXISTS `%s_ID_index` ON `%s` ( ID )" % ( args.table, args.table ))
    #db.commit()
    print( "++ Success.\n" )

# It is always good to let the user know what we are doing
# Both to confirm we're doing what they want, and because some of the steps
# can take some time (so we want to give them confidence it's still working).
# So we print out lots of messages!
print( "++ gff_to_sqlite_python_version.py" )
args = parse_arguments()
print( "++ Converting %s to %s:%s...\n" % ( args.input, args.output, args.table ))
process( args )
print( "++ Thank you for using gff_to_sqlite_python_version.py" )
