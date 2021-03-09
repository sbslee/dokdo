import os
import dokdo

def count_reads(fastq_path, delimiter='\t'):
    """Count the number of sequence reads from FASTQ.

    This command outputs two columns corresponding to the file name and
    read count, respectively.

    Parameters
    ----------
    fastq_path : str
        Path to the input FASTQ file or to the input directory
        containing FASTQ files.
    delimiter : str, default: '\t'
        Delimiter used to separate the file name and read count.
    """
    fastq_files = []

    if os.path.isdir(fastq_path):
        for r, d, f in os.walk(fastq_path):
            for fn in f:
                if '.fastq' not in fn:
                    continue
                else:
                    fastq_files.append(fn)
            break
    else:
        fastq_files.append(fastq_path)

    for fastq_file in fastq_files:
        read_count = dokdo.count_reads_one_file(fastq_file)
        print(fastq_file + delimiter + str(read_count))
