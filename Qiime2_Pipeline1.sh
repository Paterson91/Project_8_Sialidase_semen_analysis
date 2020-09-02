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
  #DADA2 contains a quality filtering step

#Clustering - clustered into Operational Taxonomic Units (OTUs)
  #Dereplication - 100% OTUs, essential first step
  #Vsearch
    #de novo -
qiime vsearch cluster-features-de-novo \
--i-table Sialidase_FeatureTable.qza \
--i-sequences Sialidase_FeatureData.qza \
--p-perc-identity 0.99 \
--o-clustered-table Sialidase_FeatureTable_de_novo.qza \
--o-clustered-sequences Sialidase_FeatureData_de_novo.qza
    #closed reference -
qiime vsearch cluster-features-closed-reference \
--i-table Sialidase_FeatureTable.qza \
--i-sequences Sialidase_FeatureData.qza \
--i-reference-sequences 85_otus.qza \
--p-perc-identity 0.85 \
--o-clustered-table Sialidase_FeatureTable_closed_0.85.qza \
--o-clustered-sequences Sialidase_FeatureData_closed_0.85.qza \
--o-unmatched-sequences Sialidase_Unmatched_closed_0.85.qza
    #open reference -
qiime vsearch cluster-features-open-reference \
--i-table Sialidase_FeatureTable.qza \
--i-sequences Sialidase_FeatureData.qza \
--i-reference-sequences 85_otus.qza \
--p-perc-identity 0.85 \
--o-clustered-table table-or-85.qza \
--o-clustered-sequences rep-seqs-or-85.qza \
--o-new-reference-sequences new-ref-seqs-or-85.qza
  #Leads to a final FeatureTable or FeatureData output dataframe

#Chimera filtering
  #DADA2 contains a Chimera filtering step



#Abundance filtering

#Downstream analysis
