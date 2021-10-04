"""
gff_to_sqlite.py
This program reads data in GFF3 format and stores it in a table in a sqlite file, indexed by the ID.
Use it like this:

    python gff_to_sqlite.py --input <path to gff file> --output <path to output file> --analysis <name>

This version has been updated to read rows in a limited number of chunks, so that it does not use excessive memory.

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

    # This version outputs two tables (gff_data and sequences)
    # and these are hard-coded for simplicoty.
    # So the --table argument is removed.
    
    parser.add_argument(
        '--overwrite',
        action = "store_true",
        help ='If specified, overwrite the table with this data.  Otherwise data will be appended.'
    )
    return parser.parse_args()

def process( args ):
    print( "++ Writing genes data from %s to %s...\n" % ( args.input, args.output ) )
    input_gff3 = open( args.input )
    output_db = sqlite3.connect( args.output )
    number_of_records = parse_gff3_by_line_to_db( input_gff3, output_db, args.overwrite )
    print( "++ ok, stored %d records in total." % number_of_records )
    input_gff3.close()
    output_db.close()
    print( "++ Success.\n" )

def parse_gff3_by_line_to_db( file, db, overwrite = False ):
    queries = build_queries( args )
    if overwrite:
        db.execute( queries['drop_gff_data'] )
    db.execute( queries['create_gff_data'] )

    chunk_size = 100000
    data = [None] * chunk_size
    index = 0
    total = 0
    for line in file:
        if line[0] == '#':
            pass
        else:
            data[index] = gff.from_gff3_line_to_dict( line )
            index = index + 1
            total = total + 1
            # write data if we have filled the chunk
            if index == chunk_size:
                db.executemany( queries['insert_gff_data'], data )
                db.commit()
                print( "  ...stored %d records..." % total )
                index = 0
                return total

    # Store the data from the last chunk, if not complete
    if index < chunk_size:
        db.executemany( queries['insert_gff_data'], data[0:index] )
        db.commit()
    return total

def build_queries( args ):
    return {
        "create_gff_data": """CREATE TABLE IF NOT EXISTS `gff_data` (
            `analysis` TEXT NOT NULL,
            `ID` TEXT,
            `Parent` TEXT,
            `Name` TEXT,
            `biotype` TEXT,
            `seqid` TEXT,
            `source` TEXT,
            `type` TEXT,
            `start` INTEGER,
            `end` INTEGER,
            `score` REAL,
            `strand` TEXT,
            `phase` TEXT,
            `attributes` TEXT
        ) ;""",
        "drop_gff_data": """DROP TABLE IF EXISTS `gff_data`""",
        "insert_gff_data": """INSERT INTO `gff_data`
VALUES( '%s', :ID, :Parent, :Name, :biotype, :seqid, :source, :type, :start, :end, :score, :strand, :phase, :attributes )
""" % ( args.analysis )
    }

# It is always good to let the user know what we are doing
# Both to confirm we're doing what they want, and because some of the steps
# can take some time (so we want to give them confidence it's still working).
# So we print out lots of messages!.print( "++ gff_to_sqlite_py" )
args = parse_arguments()
print( "++ Converting %s to %s...\n" % ( args.input, args.output ))
process( args )
print( "++ Thank you for using gff_to_sqlite.py" )
