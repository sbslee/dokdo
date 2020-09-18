import os
import zipfile
import shutil

def action_tax2seq(i_rep_seqs, i_tax, p_tax, **kwargs):
    selected = []

    with zipfile.ZipFile(i_tax, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]

    with open(zip_dir + '/data/taxonomy.tsv') as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            feature_id = fields[0]
            taxon = fields[1].split('; ')
            confidence = fields[2]

            if taxon[-1] == p_tax:
                selected.append(fields)

    fastq = {}

    with zipfile.ZipFile(i_rep_seqs, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]

    with open(zip_dir + '/data/dna-sequences.fasta') as f:
        for i, line in enumerate(f):
            if i % 2 == 0:
                feature_id = line.strip().replace('>', '')
                fastq[feature_id] = None
            else:
                seq = line.strip()
                fastq[feature_id] = seq

    for fields in selected:
        feature_id = fields[0]
        fields.append(fastq[feature_id])
        print('\t'.join(fields))
