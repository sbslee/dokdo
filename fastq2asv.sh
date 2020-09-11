#!/bin/bash

DOKDO="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

I_PATH=$1
P_TRUNC_LEN_F=$2
P_TRUNC_LEN_R=$3
P_TRIM_LEFT_F=$4
P_TRIM_LEFT_R=$5

# Make a manifest file.
python3 $DOKDO/wando.py --command make-manifest \
--i-path $I_PATH \
--o-path manifest.tsv

# Import the FASTQ files.
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path paired-end-demux.qza

# Create visualizations for the FASTQ files.
qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization demux.qzv

# Run DADA2.
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f $P_TRUNC_LEN_F \
--p-trunc-len-r $P_TRUNC_LEN_R \
--p-trim-left-f $P_TRIM_LEFT_F \
--p-trim-left-r $P_TRIM_LEFT_R \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza