# Exploring Reads in IGV


In this tutorial you are going to be inspecting the alignments you have created in a program called Integrative Genomics Viewer (IGV). IGV is a very handy piece of software which allows you to visually inspect alignments across the genome. 

First you will need to download IGV. You can find the software [here.](https://software.broadinstitute.org/software/igv/download)

Next load one of your bam files. You will be prompted to download the appropriate reference; “Pfalciparum3D7” should already be available from the provided selection.

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


