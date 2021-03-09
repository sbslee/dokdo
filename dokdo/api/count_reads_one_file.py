import gzip

def count_reads_one_file(fastq_file):
    """Count the number of sequence reads from a single FASTQ file.

    Parameters
    ----------
    fastq_file : str
        Path to the FASTQ file.

    Returns
    -------
    int
        Number of sequence reads in the FASTQ file.
    """
    n = 0
    with gzip.open(fastq_file, 'rb') as f:
        for line in f:
            n += 1
    return int(n / 4)
