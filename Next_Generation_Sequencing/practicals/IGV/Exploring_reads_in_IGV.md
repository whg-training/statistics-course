# Exploring Reads in IGV

In this tutorial you are going to be inspecting the alignments you have created in a program called Integrative Genomics Viewer (IGV). IGV is a very handy piece of software which allows you to visually inspect alignments across the genome. 

First you will need to download IGV. You can find the software [here.](https://software.broadinstitute.org/software/igv/download)

Next load one of your bam files. You will be prompted to download the appropriate reference; “Pfalciparum3D7” should already be available from the provided selection.  For this tutorial you should load (at least):

* `QG0033-C.bam`
* `QG0041-C.bam`

(You are also welcome to load the other three.)

To see reads you will need to zoom into a location on the genome.  For example, let's zoom into the gene that harbours the famous chloroquine resistance mutation, [Chloroquine resistance transporter (*CRT*)](https://plasmodb.org/plasmo/app/record/gene/PF3D7_0709000).  To search for it you will need its ID, which is `PF3D7_0709000`.  Voila!  Some reads.

IGV has lots of options that let you visualise reads in different ways.  Here are some things to try (most available by right-clicking on the reads panel).

* These reads are from paired-end sequencing.  If you right-click and choose 'view as pairs', IGV will show you this structure.
* Often there are too many reads to see.  Try using 'squished' mode to squish them up.
* On the other hand the `genes` track can be annoying when it is squished.  Right-click and choose 'expanded' to see the genes.

#### Looking at SNPs

**Question.** The mutation that causes [chloroquine resistance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2954758/) is CRT K76>T, at chr7:403,625.  (It actually involves a small haplotype including other nearby mutations).  

* Can you find these mutations?
* Parasites with the non-reference allele are resistant to chloroquine - are your parasites resistant or susceptible?
* Zoom into this location to look at the sequence.  Can you figure out the amino acid change(s)?
* How many reads support this in each sample?  What alleles do they carry?  What are their IDs?  Are they well mapped?  Are they properly paired?  Where is their pair?

**Note.** By default you are looking at all reads (including those that didn't map very well.)  There's an option hidden away under `Preferences` -> `Alignment` that lets you set a mapping quality threshold.  A good default value to put here is 20 (but remember that by doing this you are filtering some reads!)

**Note.** 'properly paired' means that the two reads in the pair align to the same region in the right orientation, and roughly fall within the distribution of insert sizes inferred from all the read pairs in the data.  The aligner (`bwa` in this case) sets a flag in the BAM file to reflect this.

#### Looking at insertions and deletions

**Question.** The variants above are single nucleotide polymorphisms (SNPs) - or rather multiple nucleotide polymorphisms (MNPs) since several come together.  What about other types of variation?  Scroll around a bit and see if you can find a deletion.  (If stuck, try looking around the gene [`PF3D7_0215300`](https://plasmodb.org/plasmo/app/search?q=PF3D7_0220300).  What read evidence supports this deletion?

**Note.** There might also be insertions, which are denoted by a purple bar across a read in IGV.  Can you find any of these?

#### The trouble with mapping is...

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

#### Examining larger structural variants

Let's have a look at a larger structural variant.  Zoom your IGV to the region 



Have a fiddle with viewing options by right-clicking on the reads and selecting from the drop-down menu. Your data was sequenced using paired end technology and by selecting “View as pairs” you’ll be shown reads joined in their respective couples. If your sample resembles the reference sequence, read pairs will fall in order along the reference with roughly a similar length of sequence between them (known as the “insert size”). However, sample genomes often don’t resemble the arbitrary reference sequences we use. Insertions, deletions and rearrangements can be inferred from read pairs that align in the “wrong” order, or with greater or smaller than expected insert sizes. This can help us define the structure of the genome in our samples. Try grouping alignments by read-pair orientation (I also like to sort by start location but try a few options and see what you like).

### **Task 1**: Have a look at your bam files and answer the following questions. 

Questions:
1.	Find a region of the genome where nothing aligns and note down the co-ordinates. Zoom into that region – why do you think no reads have aligned there? 
1.	Find a read where its pair is mapped to a different chromosome and note down its co-ordinates – can you think of reasons why this might happen?
1.  Find an example of a split read and a split read pair and note down their co-ordinates.
1.	Find an example of (i) a SNP and (ii) an indel in your reads compared to the reference and note down their co-ordinates.

We are interested in a particular region of chromosome 11 in the P. falciparum genome as there is evidence of a copy number variant (CNV) in previously published long read data.

### **Task 2**: use IGV to zoom into Pf3D7_11_v3:1,053,625-1,060,268 in your bam files. 

Compare your bam files across this region - you may notice that some samples look a bit different. See if you can answer the following questions. You may want to use [this link](https://software.broadinstitute.org/software/igv/interpreting_insert_size) and [this link](https://software.broadinstitute.org/software/igv/interpreting_pair_orientations) to help you.

Questions:
1.	Which of your alignments show evidence of a CNV in this region?
1.	In these samples do you see evidence of any deletions? How many? 
1.  Zoom into the edges of segments where you think deletions have occurred, can you tell where the breakpoints are? Do these differ between CNV samples? 
(P.S reads which only partially align to the reference are called "clipped reads" - use the internet to work out the difference between "soft" and "hard" clipping, which type do you think is being used at the edges of your deletions?)
1.	In these samples do you see any evidence of duplications? Again see if you can work out where the breakpoints are.
1.	Can you infer the structure of the CNV(s) present in this region?


