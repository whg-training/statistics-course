import gff, io

test_data = """##gff-version 3
#description: test data
chr1\tme\tgene\t1\t1000\t.\t+\t.\tID=gene1;other_data=stuff
chr1\tme\texon\t1\t1000\t.\t+\t.\tID=gene1.1;Parent=gene1
chr10\tme\tgene\t1\t1000\t.\t+\t.\tID=gene2;gene_id=my_test_gene
chr10\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript1;Parent=gene2
chr10\tme\texon\t1\t1000\t.\t+\t.\tID=my_test_exon;Parent=transcript1
chr10\tme\ttranscript\t1\t1000\t.\t+\t.\tID=transcript2;Parent=gene2
"""

X = gff.parse_gff3_to_dataframe( io.StringIO( test_data ))

assert X.loc[0]['ID'] == 'gene1'
assert X.loc[1]['ID'] == 'gene1.1'
assert X.loc[1]['Parent'] == 'gene1'
assert X.loc[2]['start'] == 1
assert X.loc[2]['end'] == 1000
assert X.loc[2]['end'] == 1000

Y = gff.parse_gff3_to_list( io.StringIO( test_data ))

assert Y[0]['ID'] == 'gene1'
assert Y[1]['ID'] == 'gene1.1'
assert Y[1]['Parent'] == 'gene1'
assert Y[2]['start'] == 1
assert Y[2]['end'] == 1000
assert Y[2]['end'] == 1000

print( "++ All tests passed! ")
