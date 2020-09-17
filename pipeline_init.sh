#!/bin/bash

X_PDP=$1
I_FASTQ=$2

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

##### Make a manifest file.
python3 $X_PDP/wando.py --command make-manifest \
--i-fastq $I_FASTQ \
--o-path manifest.tsv

echo "Saved ManifestData to: manifest.tsv"

##### Import the FASTQ files.
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path paired-end-demux.qza

##### Create visualization for the FASTQ files.
qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization demux.qzv