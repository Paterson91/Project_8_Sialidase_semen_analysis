#!/bin/bash

#Import Data

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Metadata.tsv \
  --output-path Sialidase_input.qza \
  --input-format PairedEndFastqManifestPhred33V2
#Demultiplex
  #Not required

#Cut adapters - Done prior to input

#Denoising
  #Either using DADA2 (Changes) or Deblur (Removes)
  #Reducing sequence errors and dereplicate sequences
  #By and large an arbritary choice - dada2 selected

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs Sialidase_input.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table Sialidase_FeatureTable.qza \
  --o-representative-sequences Sialidase_FeatureData.qza \
  --o-denoising-stats Sialidase_DenoisingStats.qza

#Basic quality score based filtering

#Dereplicating

qiime vsearch dereplicate-sequences \
  --i-sequences
  --o-dereplicated-table
  --o-dereplicated-sequences




#Clustering - clustered into Operational Taxonomic Units (OTUs)
  #Dereplication - 100% OTUs, essential first step
  #Vsearch
    #de novo -
    #closed reference -
    #open reference -
  #Leads to a final FeatureTable or FeatureData output dataframe








#Chimera filtering

#Abundance filtering

#Downstream analysis
