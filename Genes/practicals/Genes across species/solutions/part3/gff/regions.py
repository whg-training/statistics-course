def compute_union_of_regions( regions ):
    """Compute a non-overlapping set of regions covering the same position as the input regions.
    Input is a list of lists or tuples [ [a1,b1], ... ]
    Output is a similar list.  All regions are assumed to be closed i.e. to contain their endpoints"""
    result = []
    regions = sorted( regions, key = lambda w: w[0] )
    if len( regions ) == 0:
        return result ;
    current_region = regions[0]
    for i in range(1, len(regions)):
        current_endpoint = current_region[1]
        if regions[i][0] <= current_endpoint+1:
            current_region[1] = max( current_endpoint, regions[i][1] )
        else:
            result.append( current_region )
            current_region = regions[i]
    result.append( current_region )
    return result

def compute_genome_bases_covered( regions, sequences ):
   """Given a set of regions (as a dataframe with analysis, seqid, and start, and end columns),
   and a set of sequences (as a dataframe with analysis, seqid and sequence_length columns, return
   a dataframe showing the total number and proportion of sequence bases covered by the regions in
   each analysis."""
   import pandas

   # Utility function to add up lengths of some regions
   # Using compute_union_of_regions to avoid double-counting...
   def sum_region_lengths( regions ):
      # We convert from pandas to a list...
      aslist = regions[['start', 'end']].values.tolist()
      # ...compute the union using your function...
      union = compute_union_of_regions( aslist )
      # ...and sum to get the result.
      result = sum( [ (region[1]-region[0]+1) for region in union ] )
      # Finally, we return as a pandas.Series() object
      # because this lets us name the output variable.
      return pandas.Series(
         result,
         index = ['bases_covered']
      )
   # Utility function to sum lengths per species per chromosome
   def compute_bases_covered_in_one_chromosome( regions ):
      return (
         regions
         .groupby( [ "analysis", "seqid" ])
         .apply( sum_region_lengths )
      )
   # Utility function to add sequence lengths to the above
   def add_sequence_lengths( coverage, sequences ):
      return pandas.merge(
         coverage,
         sequences,
         left_on = [ "analysis", "seqid" ],
         right_on = [ "analysis", "seqid" ]
      )
   # Utility function to sum over chromosomes for each species...
   def sum_over_chromosomes( coverage ):
      result = coverage.groupby( "analysis" ).agg(
         bases_covered = pandas.NamedAgg( column = "bases_covered", aggfunc = sum ),
         sequence_length = pandas.NamedAgg( column = "sequence_length", aggfunc = sum )
      )
      result['proportion'] = result['bases_covered'] / result['sequence_length']
      return result
   # Now do the computation.
   per_chromosome = add_sequence_lengths(
      compute_bases_covered_in_one_chromosome( regions ),
      sequences
   )
   return sum_over_chromosomes( per_chromosome )
