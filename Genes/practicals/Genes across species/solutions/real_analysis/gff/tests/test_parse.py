import io, pandas
import unittest
import gff

class TestGff(unittest.TestCase):
    def test_parse_gff3_to_dataframe( self ): 
        test_data = """##gff-version 3
#description: test data
chr1\tme\tgene\t1\t1000\t.\t+\t.\tID=gene1;other_data=stuff
chr1\tme\texon\t1\t1000\t.\t+\t.\tID=gene1.1;Parent=gene1
chr10\tme\tgene\t1\t1000\t.\t+\t.\tID=gene2;gene_id=my_test_gene
chr10\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript1;Parent=gene2
chr10\tme\texon\t1\t1000\t.\t+\t.\tID=my_test_exon;Parent=transcript1
chr10\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript2;Parent=gene2
        """

        X = gff.parse.gff3_to_dataframe( io.StringIO( test_data ))

        self.assertEqual( X.loc[0]['ID'], 'gene1' )
        self.assertEqual( X.loc[1]['ID'], 'gene1.1' )
        self.assertEqual( X.loc[1]['Parent'], 'gene1' )
        self.assertEqual( X.loc[2]['start'], 1 )
        self.assertEqual( X.loc[2]['end'], 1000 )
        self.assertEqual( X.loc[2]['end'], 1000 )

    def test_extract_attributes_to_columns( self ):
        X = pandas.DataFrame(
            {
                "v1": [ 1, 2, 3, 4, 5 ],
                "attributes": [ "ID=1", "ID=2;Parent=1", "ID=a_string", "ID=a_string;Parent=a_string", "other=1;ID=3"]
            }
        )
        gff.parse._extract_attributes_to_columns(X, [ "ID", "Parent" ])
        self.assertEqual( X['ID'][0], "1" )
        self.assertEqual( X['ID'][1], "2" )
        self.assertEqual( X['ID'][2], "a_string" )
        self.assertEqual( X['ID'][3], "a_string" )
        self.assertEqual( X['ID'][4], "3" )
        self.assertEqual( X['Parent'][0], None )
        self.assertEqual( X['Parent'][1], "1" )
        self.assertEqual( X['Parent'][2], None )
        self.assertEqual( X['Parent'][3], "a_string" )
        self.assertEqual( X['Parent'][4], None )
        self.assertEqual( X['attributes'][0], "" )
        self.assertEqual( X['attributes'][1], "" )
        self.assertEqual( X['attributes'][2], "" )
        self.assertEqual( X['attributes'][3], "" )
        self.assertEqual( X['attributes'][4], "other=1" )

    def test_parse_attributes( self ):
        self.assertEqual( gff.parse._parse_attributes( "ID=1" )['ID'], "1" )
        self.assertEqual( gff.parse._parse_attributes( "ID=1" ).get('Parent', None), None )
        self.assertEqual( gff.parse._parse_attributes( "ID=1;Parent=5" )['ID'], "1" )
        self.assertEqual( gff.parse._parse_attributes( "ID=1;Parent=5" )['Parent'], "5" )
        self.assertEqual( gff.parse._parse_attributes( "Parent=5;ID=1" )['ID'], "1" )
