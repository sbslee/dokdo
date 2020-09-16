#!/bin/bash

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

DOKDO="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

P_CLASSIFIER=$1

##### Create the summary visualizations.
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

##### Build a phylogenetic tree.

## Carry out a multiple seqeunce alignment using Mafft.
qiime alignment mafft \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza

## Mask the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

## Create the tree using the Fasttree program.
qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza

## Root the tree using the longest root.
qiime phylogeny midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

##### Get the target sequencing depth.
DEPTH=$(python3 $DOKDO/wando.py --command compute-table-stat --i-zip table.qzv)

##### Perform the alpha rarefaction.
qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth $DEPTH \
--m-metadata-file sample-metadata.tsv \
--o-visualization alpha-rarefaction.qzv

##### Run core metrics.
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth $DEPTH \
--m-metadata-file sample-metadata.tsv \
--output-dir core-metrics-results

##### Alpha Diversity.
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization core-metrics-results/faith_pd_group-significance.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization core-metrics-results/evenness_group-significance.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/shannon_vector.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization core-metrics-results/shannon_group-significance.qzv

##### Beta diversity.
qiime emperor plot \
--i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization core-metrics-results/bray_curtis_emperor.qzv

##### Taxonomy assignment.
qiime feature-classifier classify-sklearn \
--i-classifier $P_CLASSIFIER \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization taxa-bar-plots.qzv

##### Create a Jupyter Notebook report.
python3 $DOKDO/wando.py \
--command create-report \
--o-path report.ipynb
