# Exploring Reads in IGV

## Getting started

In this tutorial you are going to be inspecting the alignments you have created in a program called Integrative Genomics Viewer (IGV). IGV is a very handy piece of software which allows you to visually inspect alignments across the genome. 

First you will need to download IGV. You can find the software [here.](https://software.broadinstitute.org/software/igv/download)

Next load one of your bam files. You will be prompted to download the appropriate reference; “Pfalciparum3D7” should already be available from the provided selection.  For this tutorial you should load (at least):

* `QG0033-C.bam`
* `QG0041-C.bam`

(You are also welcome to load the other three.)

If you followed the [snakemake tutorial](), you should also be able to load:

* Your coverage (`.bedgraph`) files for the above samples (this isn't very important because IGV shows you coverage anyway.)
* Your `variants.vcf.gz` file containing the variant calls from Octopus.

Load these in now.  Here are a few things to try:

## Looking at reads.

To see reads you will need to zoom into a location on the genome.  For example, let's zoom into the gene that harbours the famous chloroquine resistance mutation, [Chloroquine resistance transporter (*CRT*)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0709000).  To search for it you will need its ID, which is `PF3D7_0709000`.  Voila!  Some reads.

IGV has lots of options that let you visualise reads in different ways.  Here are some things to try (most available by right-clicking on the reads panel).

* These reads are from paired-end sequencing.  If you right-click and choose 'view as pairs', IGV will show you this structure.
* Often there are too many reads to see.  Try using 'squished' mode to squish them up.
* On the other hand the `genes` track can be annoying when it is squished.  Right-click and choose 'expanded' to see the genes.

## Looking at SNPs

**Question.** The mutation that causes [chloroquine resistance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2954758/) is CRT K76>T, at chr7:403,625.  (It actually involves a small haplotype including other nearby mutations).  

* Can you find these mutations?
* Parasites with the non-reference allele are resistant to chloroquine - are your parasites resistant or susceptible?
* Zoom into this location to look at the sequence.  Can you figure out the amino acid change(s)?
* How many reads support this in each sample?  What alleles do they carry?  What are their IDs?  Are they well mapped?  Are they properly paired?  Where is their pair?

**Note.** By default you are looking at all reads (including those that didn't map very well.)  There's an option hidden away under `Preferences` -> `Alignment` that lets you set a mapping quality threshold.  A good default value to put here is 20 (but remember that by doing this you are filtering some reads!)

**Note.** 'properly paired' means that the two reads in the pair align to the same region in the right orientation, and roughly fall within the distribution of insert sizes inferred from all the read pairs in the data.  The aligner (`bwa` in this case) sets a flag in the BAM file to reflect this.

## Looking at insertions and deletions

**Question.** The variants above are single nucleotide polymorphisms (SNPs) - or rather multiple nucleotide polymorphisms (MNPs) since several come together.  What about other types of variation?  Scroll around a bit and see if you can find a deletion.  (If stuck, try looking around the gene [`PF3D7_0215300`](https://plasmodb.org/plasmo/app/search?q=PF3D7_0220300).  What read evidence supports this deletion?

**Note.** There might also be insertions, which are denoted by a purple bar across a read in IGV.  Can you find any of these?

## The trouble with mapping is...

**Question.** If you scroll around a bit, after a while you'll come across regions where no reads (or very few reads) seem to align.  Can you find one?  What could cause this? 

**Hint.** Try searching for the relevant piece of genome sequence using [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).  For example, I tried this:

```
$ samtools faidx data/reference/Pf3D7_v3.fa.gz Pf3D7_07_v3:342,540-342,680
>Pf3D7_07_v3:342,540-342,680
aaaaaaaaaaaaaaagttataaaactgtataccgatattataatatatatttatataagt
atttaatacaggaatattccttgaacaaaaaaagaaatataaatataaatatatatatat
gtatatatatatatatatata
```
and then pasted the above into the [Nucleotide BLAST page](https://blast.ncbi.nlm.nih.gov/Blast.cgi).  (Another version of this can be found [on PlasmoDB](https://plasmodb.org/plasmo/app/search/transcript/UnifiedBlast).  **Note.** you may need to turn off the 'low complexity region' filter to get useful results.

## Examining larger structural variants

Let's have a look at a larger structural variant.  Zoom your IGV to the region around the gene [`Pf3D7_1127000`](https://plasmodb.org/plasmo/app/record/gene/PF3D7_1127000).  And zoom out a bit so you can see the whole region including the two flanking genes.

The samples from this tutorial were chosen because some of them contain a large structural variant in this region. 

**Challenge.** Can you find this structural variant in one or more of the samples you've loaded?

**Hints.**

* One way to find structural vairants is by looking at their effects on copy number (i.e. on read coverage.)  Do any samples have unusual copy number profiles?  (NB. this is generally easier in human sequence data because in P.falciparum the coverage varies a lot anyway due to varying AT content.)

* Your reads are paired-end reads.  If the genome has a different structure to the reference, then these might not align close together in the right orientation or insert size.  Can you find reads like this?  (Hint: the 'view as pairs', 'group by' and/or 'color by' options are very useful here.)  What is their orientation - are the read pairs the 'right' way round or the 'wrong' way round?  What does this mean?  Are these single reads or are they supported by others?

## Examining structural variant breakpoints

By a 'breakpoint' we mean the points where the segments of the sequenced genome are connected together in a different way to the reference genome.  The long, wrongly oriented / wrong insert size reads above potentially span these breakpoints.

A good way to look at this is as follows.  First zoom into one of the reads in the pair. Now right-click on a read and choose 'View mate region in split screen'.  Voila!  You are now looking at both ends of the same read pair.

To take this further there is an additional option.  If the read in the pair actually crosses a structural variant breakpoint, it could be that the end of the read does not align well.  If so it will typically be 'soft-clipped' by the aligner.  (The [`CIGAR string`](https://sites.google.com/site/bioinformaticsremarks/bioinfo/sam-bam-format/what-is-a-cigar) for the alignment will have a sequence of 'S's at the end).  Under `Preferences` -> `Alignments` there is an option 'Show soft-clipped bases' to visualise these.  Go and turn this on.  Does this help figure out the basepair location of the possible breakpoint?

**NB.** In the current data this works best for the reads at `Pf3D7_11_v3:1,058,745-1,058,915` - if you processed the [full coverage data]() then it might work at all the breakpoints.

## Figuring out the structural variant structure

**Extra hard challenge.** Can you figure out the structure of the structural variant genomes from these reads?

**Warning.** This really is hard!  If you want some inspiration, see Figure 6 of [this paper on the Dantu blood group](https://doi.org/10.1126/science.aam6393) to see a similar process applied to a human structural variant.

## Summary

Congratulations!  You are now an IGV expert!
