#!/bin/bash

P_TRUNC_LEN_F=$1
P_TRUNC_LEN_R=$2
P_TRIM_LEFT_F=$3
P_TRIM_LEFT_R=$4
P_N_THREADS=$5

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

##### Run DADA2.
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f $P_TRUNC_LEN_F \
--p-trunc-len-r $P_TRUNC_LEN_R \
--p-trim-left-f $P_TRIM_LEFT_F \
--p-trim-left-r $P_TRIM_LEFT_R \
--p-n-threads $P_N_THREADS \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

##### Create the visualization files.
# qiime metadata tabulate \
# --m-input-file denoising-stats.qza \
# --o-visualization denoising-stats.qzv