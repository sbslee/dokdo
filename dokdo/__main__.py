import argparse
from .cli import commands
from .version import __version__

def main():
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version and exit."
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        help="Name of the command."
    )

    subparsers.required = True

    collapse_parser = subparsers.add_parser(
        "collapse",
        add_help=False,
        description=("This command creates seven collapsed ASV tables, one "
                     "for each taxonomic level (i.e. `level-1.csv` to "
                     "`level-7.csv`)."),
        help=("This command creates seven collapsed ASV tables, one for each "
              "taxonomic level."),
    )

    collapse_parser._optionals.title = "Arguments"

    collapse_parser.add_argument(
        "-t",
        "--table",
        metavar="PATH",
        required=True,
        help="Path to the input table.qza file. [required]"
    )
    collapse_parser.add_argument(
        "-x",
        "--taxonomy",
        metavar="PATH",
        required=True,
        help="Path to the input taxonomy.qza file. [required]"
    )
    collapse_parser.add_argument(
        "-o",
        "--output-dir",
        metavar="PATH",
        required=True,
        help=("Path to the output directory. [required]")
    )

    collapse_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )







    make_manifest_parser = subparsers.add_parser(
        'make-manifest',
        description=("This command creates a manifest file from a directory "
                     "containing FASTQ files. The file names must include "
                     "either '_R1_001.fastq' or '_R1_002.fastq'. The word "
                     "before the third-to-last underscore will be set as the "
                     "sample ID. For example, a file named "
                     "'EXAMPLE_S1_R1_001.fastq.gz' will produce 'EXAMPLE' as "
                     "the sample ID and 'EXAM_PLE_S1_R1_001.fastq.gz', "
                     "'EXAM_PLE'."),

        help=("This command creates a manifest file from a directory "
              "containing FASTQ files."),
    )
    make_manifest_parser.add_argument(
        'fastq_dir',
        help="Path to the directory containing input FASTQ files (.fastq.gz)."
    )
    make_manifest_parser.add_argument(
        'output',
        help="Path to the output manifest file (.tsv)."
    )










    add_metadata_parser = subparsers.add_parser(
        'add-metadata',
        description=("This command adds new columns to an existing "
                     "sample-metadata file. The 'metadata' file and the "
                     "'columns' file must have at least one overlapping "
                     "column."),
        help=("This command adds new columns to an existing sample-metadata "
              "file."),
    )
    add_metadata_parser.add_argument(
        'metadata',
        help="Path to the input sample-metadata.tsv file."
    )
    add_metadata_parser.add_argument(
        'columns',
        help=("Path to a file containing the new columns to be added (.tsv)."
              "The first row should be column names.")
    )
    add_metadata_parser.add_argument(
        'output',
        help="Path to the output sample-metadata.tsv file."
    )










    merge_metadata_parser = subparsers.add_parser(
        'merge-metadata',
        description=("This command merges two or more sample-metadata.tsv "
                     "files vertically. All files are assumed to have the "
                     "same column names."),
        help=("This command merges two or more sample-metadata.tsv files "
              "vertically."),
    )
    merge_metadata_parser.add_argument(
        'metadata',
        nargs='+',
        help="Paths to the sample-metadata.tsv files to be merged.",
    )
    merge_metadata_parser.add_argument(
        'output',
        help="Path to the output sample-metadata.tsv file.",
    )










    summarize_parser = subparsers.add_parser(
        'summarize',
        description=("This command extracts summary statistics from an "
                     "Artifact file with the semantic type "
                     "`FeatureTable[Frequency]`."),
        help=("This command extracts summary statistics for an Artifact file."),
    )
    summarize_parser.add_argument(
        'input',
        help="Path to the input Artifact file.",
    )
    summarize_parser.add_argument(
        '-v', '--verbose', action='store_true',
        help="Print a verbose version of the results.",
    )









    args = parser.parse_args()
    command = args.command
    delattr(args, "command")
    commands[command](**vars(args))

if __name__ == "__main__":
    main()
