import io, pandas, unittest
import gff

class TestGff(unittest.TestCase):
    def test_count_exons_per_transcript( self ):
        # Set up some test data
        transcript_data = """##gff-version 3
#description: test data
chr1\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript1;Parent=gene1
chr1\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript2;Parent=gene1
"""
        exon_data = """##gff-version 3
#description: test data
chr1\tme\texon\t1\t400\t.\t+\t.\tID=exon1;Parent=transcript1
chr1\tme\texon\t500\t1000\t.\t+\t.\tID=exon2;Parent=transcript1
"""
        transcripts = gff.parse.gff3_to_dataframe( io.StringIO( transcript_data ))
        exons = gff.parse.gff3_to_dataframe( io.StringIO( exon_data ))
        summary = gff.summary.count_exons_per_transcript( transcripts, exons )
        self.assertEqual( summary['ID'][0], 'transcript1' )
        self.assertEqual( summary['ID'][1], 'transcript2' )
        self.assertEqual( summary['Parent'][0], 'gene1' )
        self.assertEqual( summary['Parent'][1], 'gene1' )
        self.assertEqual( summary['number_of_exons'][0], 2 )
        self.assertEqual( summary['number_of_exons'][1], 0 )

    def test_summarise_genes( self ):
        # Set up some test data
        gene_data = """##gff-version 3
#description: test data
chr1\tme\tgene\t1\t1000\t.\t+\t.\tID=gene1
chr1\tme\tgene\t1\t1000\t.\t+\t.\tID=gene2
"""
        transcript_data = """##gff-version 3
#description: test data
chr1\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript1;Parent=gene1
chr1\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript2;Parent=gene1
"""
        exon_data = """##gff-version 3
#description: test data
chr1\tme\texon\t1\t400\t.\t+\t.\tID=exon1;Parent=transcript1
chr1\tme\texon\t500\t1000\t.\t+\t.\tID=exon2;Parent=transcript1
"""
        genes = gff.parse.gff3_to_dataframe( io.StringIO( gene_data ))
        transcripts = gff.parse.gff3_to_dataframe( io.StringIO( transcript_data ))
        exons = gff.parse.gff3_to_dataframe( io.StringIO( exon_data ))

        summary = gff.summary.summarise_genes( genes, transcripts, exons )
        # Or use the python version:
        # summary = gff.summarise_genes_python_version( genes, transcripts, exons )
        self.assertEqual( summary['ID'][0], 'gene1' )
        self.assertEqual( summary['ID'][1], 'gene2' )
        self.assertEqual( summary['number_of_transcripts'][0], 2 )
        self.assertEqual( summary['number_of_transcripts'][1], 0 )
        self.assertEqual( summary['average_number_of_exons'][0], 1 )
        import math
        self.assertTrue( math.isnan( summary['average_number_of_exons'][1] ))

if __name__ == '__main__':
    unittest.main()
