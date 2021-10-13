## A list of relevant GFF files.

A list of gff files from ensembl can be found in [`solutions/ensembl_gff_listing.txt`](solutions/ensembl_gff_listing.txt).

If you look [on ensembl](http://ftp.ensembl.org/pub/current_gff3/) you will find lots of these gff files - I have selected some as follows:

* if a species has a file of the form `[species].[version].chr.gff3.gz` then I used that.  
* otherwise, I used one of the form `[species].[version].gff3.gz`.

The idea is to get annotations for the whole genome.  Genome assemblies typically have a set of core 'primary' contigs, and then might have a set of additional contigs.
I *think* (but not sure) that the `.chr.gff.gz` are the annotations for the main reference chromosomes, so have picked these.

You should read the README files for the annotations for details.
