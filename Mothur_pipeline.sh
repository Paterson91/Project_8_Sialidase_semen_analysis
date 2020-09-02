#Mothur Commands - Sialidase Exp

#Input files (.fastq format)
make.file(inputdir=., type=fastq, prefix=stability)

#Combining two sets of reads  for each sample
make.contigs(file=stability.files, processors=100)

#Resulting report
summary.seqs(fasta=stability.trim.contigs.fasta)

#Screen sequences - QC Check
  #Note - Sequences should be less that Xbp
  #MiSeq max length should be 2x150bp
screen.seqs(fasta=stability.trim.contigs.fasta,group=stability.contigs.groups,maxambig=0,maxlength=465)

#Select unique sequences
unique.seqs(fasta=stability.trim.contigs.good.fasta)

#Make count of uniquely selected sequences
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)

#Trim reference sequences to adapter regions
trim.seqs(fasta=/home/ap14958/FASTA/SILVA/SILVA_132_SSURef_Nr99_tax_silva.fasta, oligos=oligos.oligo)


#Align to reference sequences (Selected for V4 region)
#Forward Primer; CCTACGGGNGGCWGCAG
#Reverse Primer; GACTACHVGGGTATCTAATCC (Reverse Complement; GGATTAGATACCCGTAGTC)
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=100)
pcr.seqs(fasta=stability.trim.contigs.good.unique.fasta,processors=100)
