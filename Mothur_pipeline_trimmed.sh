#Mothur Commands - Sialidase Exp

#Input files (.fastq format)
make.file(inputdir=., type=fastq, prefix=stability)

#Combining two sets of reads for each sample
make.contigs(file=stability.files, processors=100)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.fasta
#/mnt/Scratch/Sialidase/Trimmed/stability.scrap.contigs.fasta
#/mnt/Scratch/Sialidase/Trimmed/stability.contigs.report
#/mnt/Scratch/Sialidase/Trimmed/stability.contigs.groups

#Resulting report
summary.seqs(fasta=stability.trim.contigs.fasta)
#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	40	40	0	3	1
#2.5%-tile:	1	285	285	0	4	240605
#25%-tile:	1	286	286	0	5	2406041
#Median: 	1	289	289	0	6	4812081
#75%-tile:	1	313	313	0	6	7218121
#97.5%-tile:	1	465	465	11	7	9383557
#Maximum:	1	602	602	123	301	9624161
#Mean:	1	317	317	1	5
# of Seqs:	9624161
#It took 66 secs to summarize 9624161 sequences.
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.summary

#Screen sequences - QC Check
#V3-4;  451bp
screen.seqs(fasta=stability.trim.contigs.fasta,group=stability.contigs.groups,maxambig=0,minlength=400,maxlength=500)
screen.seqs(fasta=stability.trim.contigs.fasta,group=stability.contigs.groups,maxambig=0,minlength=400,maxlength=500)

#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.fasta
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.bad.accnos
#/mnt/Scratch/Sialidase/Trimmed/stability.contigs.good.groups

#Select unique sequences
unique.seqs(fasta=stability.trim.contigs.good.fasta)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.names
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.fasta

#Make count of uniquely selected sequences
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.count_table

#Download 16s rRNA sequence from NCBI
>NC_003888.3:4530650-4532177 Streptomyces coelicolor A3(2) chromosome, complete genome
CCCGATTACGGGTATACATTCACGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAAC
ACATGCAAGTCGAACGATGAACCACTTCGGTGGGGATTAGTGGCGAACGGGTGAGTAACACGTGGGCAAT
CTGCCCTTCACTCTGGGACAAGCCCTGGAAACGGGGTCTAATACCGGATACTGACCCTCGCAGGCATCTG
CGAGGTTCGAAAGCTCCGGCGGTGAAGGATGAGCCCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTC
ACCAAGGCGACGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACACTGGGACTGAGACACGGCCCAGAC
TCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGAGG
GATGACGGCCTTCGGGTTGTAAACCTCTTTCAGCAGGGAAGAAGCGAAAGTGACGGTACCTGCAGAAGAA
GCGCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTGTCCGGAATTATTGGGC
GTAAAGAGCTCGTAGGCGGCTTGTCACGTCGGTTGTGAAAGCCCGGGGCTTAACCCCGGGTCTGCAGTCG
ATACGGGCAGGCTAGAGTTCGGTAGGGGAGATCGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCA
GGAGGAACACCGGTGGCGAAGGCGGATCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGC
GAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGGCACTAGGTGTGGGCAACATTCCACGT
TGTCCGTGCCGCAGCTAACGCATTAAGTGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGG
AATTGACGGGGGCCCGCACAAGCGGCGGAGCATGTGGCTTAATTCGACGCAACGCGAAGAACCTTACCAA
GGCTTGACATACACCGGAAAGCATCAGAGATGGTGCCCCCCTTGTGGTCGGTGTACAGGTGGTGCATGGC
TGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCCGTGTTGCC
AGCAAGCCCTTCGGGGTGTTGGGGACTCACGGGAGACCGCCGGGGTCAACTCGGAGGAAGGTGGGGACGA
CGTCAAGTCATCATGCCCCTTATGTCTTGGGCTGCACACGTGCTACAATGGCCGGTACAATGAGCTGCGA
TACCGCAAGGTGGAGCGAATCTCAAAAAGCCGGTCTCAGTTCGGATTGGGGTCTGCAACTCGACCCCATG
AAGTCGGAGTCGCTAGTAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACC
GCCCGTCACGTCACGAAAGTCGGTAACACCCGAAGCCGGTGGCCCAACCCCTTGTGGGAGGGAGCTGTCG
AAGGTGGGACTGGCGATTGGGACGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGCG
#Forward Primer; Bakt_341F CCTACGGGNGGCWGCAG (regex; CCTACGGG.GGC.GCAG)
#Reverse Primer; Bakt_805R GACTACHVGGGTATCTAATCC (Reverse Complement; GGATTAGATACCC..GTAGTC)

#Create pcrTest.oligos file - BASH
printf "forward\tCCTACGGGNGGCWGCAG\nreverse\tGACTACHVGGGTATCTAATCC\n" > pcrTest.oligos

pcr.seqs(fasta=S.coelicolor.fasta, oligos=pcrTest.oligos)
#Produces the trimmed V3-4 region using the above primers (451bp)
#Did not include primers, so shorter and possibly leading to reduced alignment downstream
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/S.coelicolor.pcr.fasta

#Align to reference sequences (Selected for V3-4 region approx. 450bp)
align.seqs(fasta=S.coelicolor_v3-4.fasta, reference=silva.nr_v132.align)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/S.coelicolor_v3-4.align
#/mnt/Scratch/Sialidase/Trimmed/S.coelicolor_v3-4.align.report

summary.seqs(fasta=S.coelicolor_v3-4.align)
#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	6388	25316	445	0	4	1
#2.5%-tile:	6388	25316	445	0	4	1
#25%-tile:	6388	25316	445	0	4	1
#Median: 	6388	25316	445	0	4	1
#75%-tile:	6388	25316	445	0	4	1
#97.5%-tile:	6388	25316	445	0	4	1
#Maximum:	6388	25316	445	0	4	1
#Mean:	6388	25316	445	0	4
# of Seqs:	1
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/S.coelicolor_v3-4.summary

pcr.seqs(fasta=silva.nr_v132.align, start=6388, end=25316, keepdots=F, processors=100)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/silva.nr_v132.pcr.align

#Rename to something more useful
rename.file(input=silva.nr_v132.pcr.align, new=silva.v34.fasta)

summary.seqs(fasta=silva.v34.fasta)
#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	15382	346	0	3	1
#2.5%-tile:	1	18928	435	0	4	5328
#25%-tile:	1	18928	440	0	4	53280
#Median: 	1	18928	460	0	5	106560
#75%-tile:	1	18928	464	0	6	159840
#97.5%-tile:	1	18928	605	1	7	207792
#Maximum:	3865	18928	1666	5	16	213119
#Mean:	1	18927	467	0	5
# of Seqs:	213119
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/silva.v34.summary

align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v34.fasta, flip=t)
#It took 3399 secs to align 302921 sequences.
#[WARNING]: 18770 of your sequences generated alignments that eliminated too many bases, a list is provided in /mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.flip.accnos.
#[NOTE]: 14948 of your sequences were reversed to produce a better alignment.
#It took 3399 seconds to align 302921 sequences.
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.align
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.align.report
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.flip.accnos

summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)
#             Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	    0	     0	0	0	1	1
#2.5%-tile:	  1	     25	10	0	3	12287
#25%-tile:    1	18928	439	0	4	122864
#Median: 	    1	18928	454	0	5	245728
#75%-tile:	  1	18928	463	0	6	368591
#97.5%-tile:	46	18928	464	0	7	479168
#Maximum:	    18928	18928	484	0	15	491454
#Mean:	      255	18042	425	0	4
# of unique seqs:	302921
#total # of seqs:	491454
#It took 23 secs to summarize 491454 sequences.
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.summary

screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1, end=18928, maxhomop=8)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.summary
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.align
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.bad.accnos
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.good.count_table

summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)
#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	18928	407	0	3	1
#2.5%-tile:	1	18928	437	0	4	11334
#25%-tile:	1	18928	439	0	4	113333
#Median: 	1	18928	459	0	5	226666
#75%-tile:	1	18928	464	0	6	339999
#97.5%-tile:	1	18928	464	0	7	441998
#Maximum:	1	18928	484	0	8	453331
#Mean:	1	18928	451	0	5
# of unique seqs:	276135
#total # of seqs:	453331
#It took 20 secs to summarize 453331 sequences.
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.summary

filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
#Length of filtered alignment: 1169
#Number of columns removed: 17759
#Length of the original alignment: 18928
#Number of sequences used to construct filter: 276135
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.filter
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.filter.fasta

#Repeat unique selection to reduce redundancy
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)
#Output File Names:
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.filter.count_table
#/mnt/Scratch/Sialidase/Trimmed/stability.trim.contigs.good.unique.good.filter.unique.fasta

#Further de-noise by pre-clustering allowing up to 2 differences between sequences
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Now left with unique sequences with as little sequencing error as possible. Now to remove chimeras
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table
#stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras
#stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos

#Removal of chimeric sequences
remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
#Removed 85503 sequences from your fasta file.
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta

#What are we left with?
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	1169	411	0	3	1
#2.5%-tile:	1	1169	438	0	4	8917
#25%-tile:	1	1169	439	0	4	89166
#Median: 	1	1169	459	0	5	178331
#75%-tile:	1	1169	464	0	6	267496
#97.5%-tile:	1	1169	464	0	7	347745
#Maximum:	1	1169	467	0	8	356661
#Mean:	1	1169	452	0	5
# of unique seqs:	111663
#total # of seqs:	356661
#It took 9 secs to summarize 356661 sequences.
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary

#This will classify each sequence into a taxonomy
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=silva.nr_v132.align, taxonomy=silva.nr_v132.tax, cutoff=80)
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.tax.summary

#In order to remove undesirables that are non-bacterial
remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.accnos
#stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta

summary.tax(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)




#Assessing error rates based on mock community
#This can only be conducted if you have co-sequenced a mock community of a known composition
#Here we have no mock community of known composition





#Now assess OTUs (Operational Taxonomic Unit)
#We can now cluster sequences into OTUs
dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#You did not set a cutoff, using 0.03.
#Clustering stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist
#iter	time	label	num_otus	cutoff	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score

#0.03
#0	0	0.03	111313	0.03	0	6.09466e+09	0	1.00575e+08	0	1	0	0.983766	1	0.983766	0	0
#1	172	0.03	8724	0.03	9.72363e+07	6.09069e+09	3.97214e+06	3.33889e+06	0.966802	0.999348	0.960753	0.999452	0.960753	0.99882	0.963173	0.963768
#2	172	0.03	8220	0.03	9.73278e+07	6.09064e+09	4.02371e+06	3.2474e+06	0.967712	0.99934	0.960299	0.999467	0.960299	0.998826	0.963402	0.963991
#3	164	0.03	8210	0.03	9.73265e+07	6.09064e+09	4.02109e+06	3.24868e+06	0.967699	0.99934	0.960324	0.999467	0.960324	0.998827	0.963408	0.963997

#It took 1369 seconds to cluster
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sensspec

#Next we want to know how many sequences are in each OTU from each group and we can do this using the Make.shared command. Here we tell Mothur that we’re really only interested in the 0.03 cutoff level:
make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared


#We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the Classify.otu command:
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label=0.03)
#Output File Names:
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
#stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary

#OTU-Based Analysis
#Count how many sequences we have in each sample
count.groups(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)

BES004 contains 2636.
BES005 contains 9419.
BES006 contains 10848.
BES008 contains 3022.
BES009 contains 6273.
BES010 contains 2024.
BES011 contains 12604.
BES013 contains 12727.
BES014 contains 13272.
BES015 contains 8193.
BES016 contains 13191.
BES017 contains 9782.
BES019 contains 7886.
BES023 contains 3081.
BES024 contains 4914.
BES027 contains 7405.
BES028 contains 19749.
BES030 contains 688.
BES034 contains 16025.
BES035 contains 2075.
BES036 contains 464.
BES037 contains 662.
BES038 contains 122.
BES039 contains 530.
BES043 contains 491.
BES045 contains 331.
BES046 contains 198.
BES047 contains 371.
BES050 contains 14446.
BES057 contains 4554.
BES058 contains 6998.
BES060 contains 199.
BES061 contains 350.
BES063 contains 1482.
BES069 contains 5485.
BES070 contains 168.
BES071 contains 136.
BES072 contains 2844.
BES074 contains 761.
BES081 contains 946.
BES082 contains 27.
BES083 contains 4196.
BES084 contains 129.
BES093 contains 8478.
BES097 contains 8428.
BES100 contains 5959.
BES101 contains 2689.
BES102 contains 2291.
BES105 contains 17939.
BES107 contains 2699.
BES141 contains 4507.
BES142 contains 628.
BES143 contains 306.
BES144 contains 30.
BES145 contains 1381.
BES146 contains 6865.
BES147 contains 6003.
BES148 contains 534.
BES150 contains 1845.
BES152 contains 115.
BES153 contains 74.
BES155 contains 169.
BES157 contains 952.
BES158 contains 1763.
BES159 contains 6023.
BES160 contains 386.
BES161 contains 10728.
BES162 contains 82.
BPS01A contains 71.
BPS02B contains 194.
BPS03C contains 3192.
BPS04D contains 774.
BPS05E contains 977.
BPS06F contains 1406.
BPS07G contains 1077.
BPS08H contains 1647.
BPS10J contains 2987.
BPS11K contains 212.
BPS12L contains 1032.
BPS13M contains 77.
BPS14N contains 1794.
BPS15O contains 130.
BPS16P contains 1985.
BPS19S contains 4154.
BPS20T contains 7623.
BPS21U contains 11970.
BPS22V contains 975.
BPS24X contains 2784.
BPS26Z contains 43.
BPS27A contains 399.
BPS28B contains 871.
BPS29C contains 1321.
BPS30D contains 73.
BPS32F contains 3897.
BPS33G contains 1518.
BPS34H contains 229.

Size of smallest group: 27.

Total seqs: 355990.
