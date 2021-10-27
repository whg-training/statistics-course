import pandas
from .parse import gff3_to_dataframe
from .parse import extract_attributes_to_columns
from .sequences import sequences_from_gff3_metadata
from .summary import summarise_transcripts_and_exons_per_gene
from .regions import compute_genome_bases_covered

class UnpackedGFF:
    def __init__( self, analysis ):
        self.analysis = analysis
        self.genes = None
        self.transcripts = None
        self.exons = None
        self.cds = None
        self.sequences = None

    def summary( self ):
        return "++ UnpackedGFF.summarise(): %d genes, %d transcripts, %d exons and %d coding sequence entries; %d sequences." % (
           self.genes.shape[0],
           self.transcripts.shape[0],
           self.exons.shape[0],
           self.cds.shape[0],
           self.sequences.shape[0]
        )

    @classmethod
    def from_gff3( cls, filename, analysis ):
        result = UnpackedGFF( analysis )
        data = UnpackedGFF.load_and_sanitise_data( filename )
        
        print( "++ UnpackedGFF.from_gff3(): Extracting regular protein-coding genes..." )
        result.genes = data[ data['type'].isin( [ 'gene', 'protein_coding_gene' ] )].copy()
        extract_attributes_to_columns( result.genes, [ 'biotype' ], 3 )
        result.genes.loc[result.genes.type == 'protein_coding_gene', 'biotype' ] = 'protein_coding'
        result.genes.loc[result.genes.type == 'protein_coding_gene', 'type' ] = 'gene'
        result.genes = result.genes[ result.genes['biotype'] == 'protein_coding' ]
        result.genes.insert( 0, 'analysis', analysis )
        
        print( "++ UnpackedGFF.from_gff3(): Extracting transcripts..." )
        result.transcripts = data[ (data['type'] == 'mRNA') & (data['Parent'].isin( result.genes['ID'] ))].copy()
        result.transcripts.insert( 0, 'analysis', analysis )
        extract_attributes_to_columns( result.transcripts, [ 'tag', 'transcript_support_level' ] )
        
        print( "++ UnpackedGFF.from_gff3(): Extracting exons..." )
        result.exons = data[ (data['type'] == 'exon') & (data['Parent'].isin( result.transcripts['ID'] ))].copy()
        result.exons.insert( 0, 'analysis', analysis )
        print( "++ UnpackedGFF.from_gff3(): Extracting cds..." )
        result.cds = data[ (data.type == 'CDS') & (data['Parent'].isin( result.transcripts['ID'] ))].copy()
        result.cds.insert( 0, 'analysis', analysis )
        
        result.sequences = UnpackedGFF.load_sequences( filename )
        result.sequences.insert( 0, 'analysis', analysis )
        result.sequences.insert( result.sequences.shape[1], 'sequence_length', result.sequences['end'] - result.sequences['start'] + 1 )
        
        result.genes = result.genes[[
            "analysis",
            "ID",
            "Parent",
            "Name",
            "biotype",
            "seqid",
            "source",
            "start",
            "end",
            "strand",
            "attributes"
        ]]
        result.transcripts = result.transcripts[[
            "analysis",
            "ID",
            "Parent",
            "seqid",
            "source",
            "start",
            "end",
            "strand",
            "tag",
            "transcript_support_level",
            "attributes"
        ]]
        result.exons = result.exons[[
            "analysis",
            "ID",
            "Parent",
            "seqid",
            "source",
            "start",
            "end",
            "strand"
        ]]
        result.cds = result.cds[[
            "analysis",
            "ID",
            "Parent",
            "seqid",
            "source",
            "start",
            "end",
            "strand"
        ]]
        return result

    @classmethod
    def load_and_sanitise_data( self, filename ):
        print( "++ UnpackedGFF: Loading gene records from %s..." % filename )
        if filename.endswith( ".gz" ):
            import gzip
            fileinput = gzip.open( filename, 'r' )
        else:
            fileinput = open( filename, 'rt' )
        result = gff3_to_dataframe( fileinput )
        # Remove unwanted parts of IDs: the "gene:" and "transcript:" prefixes
        def fix_id( id ):
            if id is None:
                return None
            else:
                return id.replace( "gene:", "" ).replace( "transcript:", "" ).replace( "CDS:", "" )
        result.loc[ :, 'ID'] = result['ID'].apply( fix_id )
        result.loc[ :, 'Parent'] = result['Parent'].apply( fix_id )
        print( "++ ok, %d records loaded." % result.shape[0] )
        return result

    @classmethod
    def load_sequences( self, filename ):
        print( "++ UnpackedGFF: Loading sequences from %s..." % filename )
        if filename.endswith( ".gz" ):
            import gzip
            fileinput = gzip.open( filename, 'rt' )
        else:
            fileinput = open( filename, 'r' )
        return sequences_from_gff3_metadata( fileinput )

def summarise_genes( data ):
    import pandas
    result = {
        "per_gene": summarise_transcripts_and_exons_per_gene( data.genes, data.transcripts, data.exons ),
        "statistics": None
    }
    print( result['per_gene'] )
    result['statistics'] = result['per_gene'].groupby( 'analysis' ).agg(
        total_genes = pandas.NamedAgg(
            column = 'ID',
            aggfunc = lambda x: x.notnull().sum()
        ),
        total_single_exon_genes = pandas.NamedAgg(
            column = 'mean_exons',
            aggfunc = lambda x: ( x == 1 ).sum()
        ),
        proportion_single_exon_genes = pandas.NamedAgg(
            column = 'mean_exons',
            aggfunc = lambda x: (x == 1).sum() / x.notnull().sum()
        ),
        longest_gene = pandas.NamedAgg(
            column = "max_transcript_length",
            aggfunc = lambda x: x.max()
        ),
        shortest_gene = pandas.NamedAgg(
            column = "max_transcript_length",
            aggfunc = lambda x: x.min()
        ),
        highest_exon_count = pandas.NamedAgg(
            column = "mean_exons",
            aggfunc = lambda x: x.max()
        ),
        highest_transcript_count = pandas.NamedAgg(
            column = "number_of_transcripts",
            aggfunc = lambda x: x.max()
        )
    )
    result[ 'statistics' ].reset_index( inplace = True )
    return result

def compute_coverage_statistics( data ):
    """Compute the protein-coding gene, exon, and cds genome coverage"""
    import pandas
    
    def compute_coverages( genes, exons, cds, sequences ):
        # the .reset_index() is used here to turn the analysis column back into a normal column.
        # (otherwise it is an 'index' and behaves differently)
        result = {
            "genes": compute_genome_bases_covered( genes, sequences ).reset_index(),
            "exons": compute_genome_bases_covered( exons, sequences ).reset_index(),
            "cds": compute_genome_bases_covered( cds, sequences ).reset_index()
        }
        return result
    
    def build_single_table( coverages ):
        # Now build a single table
        result = coverages['genes'][['analysis', 'sequence_length']]
        for what in [ 'genes', 'exons', 'cds' ]:
            result = pandas.merge(
                result,
                coverages[what][['analysis', 'bases_covered', 'proportion' ]],
                left_on = 'analysis',
                right_on = 'analysis'
            )
            result.rename(
                columns = {
                    "bases_covered": "%s:bases_covered" % what,
                    "proportion": "%s:proportion_covered" % what
                },
                inplace = True
            )
        return result

    coverages = compute_coverages( data.genes, data.exons, data.cds, data.sequences )
    return build_single_table( coverages )
