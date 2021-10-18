import os
from pathlib import Path

def make_manifest(fastq_dir, output_file):
    """
    Create a manifest file (.tsv) from a directory containing FASTQ files.

    This command assumes that FASTQ filenames end with a suffix such as
    '_S0_R1_001.fastq' or '_S14_R2_001.fastq'. The word before the
    third-to-last underscore ('_') will be used as sample ID. For example, a
    file named 'EXAMPLE_S1_R1_001.fastq.gz' will set 'EXAMPLE' as sample ID.
    Undertermined reads (e.g. 'Undetermined_S0_R1_001.fastq') will not be
    included in the output file.

    Parameters
    ----------
    fastq_dir : str
        Directory containing input FASTQ files.
    output_file : str
        Manifest file.
    """
    fastq_dir = Path(fastq_dir).resolve()

    files = {}

    for r, d, f in os.walk(fastq_dir):
        for x in f:
            name = '_'.join(x.split('_')[:-3])

            if 'Undetermined' in x:
                continue

            if '_R1_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][0] = f"{r}/{x}"
            elif '_R2_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][1] = f"{r}/{x}"
            else:
                pass

    with open(output_file, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath',
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for name in sorted(files):
            fields = [name, files[name][0], files[name][1]]
            f.write('\t'.join(fields) + '\n')
