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
        "--table-file",
        metavar="PATH",
        required=True,
        help="Path to the input table.qza file. [required]"
    )

    collapse_parser.add_argument(
        "-x",
        "--taxonomy-file",
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
        "make-manifest",
        add_help=False,
        description=("This command creates a manifest file (.tsv) from a "
                     "directory containing FASTQ files. The file names must "
                     "include either '_R1_001.fastq' or '_R1_002.fastq'. "
                     "The word before the third-to-last underscore will "
                     "be set as the sample ID. For example, a file named "
                     "'EXAMPLE_S1_R1_001.fastq.gz' will produce 'EXAMPLE' as "
                     "sample ID and 'EXAM_PLE_S1_R1_001.fastq.gz', "
                     "'EXAM_PLE'."),

        help=("This command creates a manifest file (.tsv) from a directory "
              "containing FASTQ files."),
    )

    make_manifest_parser._optionals.title = "Arguments"

    make_manifest_parser.add_argument(
        "-i",
        "--fastq-dir",
        required=True,
        metavar="PATH",
        help="Path to the directory containing input FASTQ files. [required]"
    )

    make_manifest_parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        metavar="PATH",
        help="Path to the output file. [required]"
    )

    make_manifest_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
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





    summarize_parser = subparsers.add_parser(
        'summarize',
        add_help=False,
        help=("Extract summary or verbose data from an Artifact file."),
        description=("Extract summary or verbose data from an Artifact "
                     "file. This command automatically detects the input "
                     "file's semantic type and then extracts summary or "
                     "verbose data from it.")
    )

    summarize_parser._optionals.title = "Arguments"

    summarize_parser.add_argument(
        "-i",
        "--input-file",
        help="Path to the input Artifact file. [required]",
    )

    summarize_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print a verbose version of the results.",
    )

    summarize_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    prepare_lefse_parser = subparsers.add_parser(
        "prepare-lefse",
        add_help=False,
        help=("Create a text file which can be used as input for "
              "the LEfSe tool."),
        description=("Create a text file which can be used as input for the "
                     "LEfSe tool. This command 1) collapses the input "
                     "feature table at the genus level, 2) computes "
                     "relative frequency of the features, 3) performs "
                     "sample filtration if requested, 4) changes the "
                     "format of feature names, 5) adds the relevant "
                     "metadata as 'Class', 'Subclass', and 'Subject', "
                     "and 6) writes a text file which can be used as "
                     "input for LEfSe.")
    )

    prepare_lefse_parser._optionals.title = "Arguments"

    prepare_lefse_parser.add_argument(
        "-t",
        "--table-file",
        metavar="PATH",
        required=True,
        help=("Path to the table file with the 'FeatureTable[Frequency]' "
              "type. [required]")
    )

    prepare_lefse_parser.add_argument(
        "-x",
        "--taxonomy-file",
        metavar="PATH",
        required=True,
        help=("Path to the taxonomy file with the 'FeatureData[Taxonomy]' "
              "type. [required]")
    )

    prepare_lefse_parser.add_argument(
        "-m",
        "--metadata-file",
        metavar="PATH",
        required=True,
        help="Path to the metadata file. [required]"
    )

    prepare_lefse_parser.add_argument(
        "-o",
        "--output-file",
        metavar="PATH",
        required=True,
        help=("Path to the output file. [required]")
    )

    prepare_lefse_parser.add_argument(
        "-c",
        "--class-col",
        metavar="TEXT",
        required=True,
        help="Metadata column used as 'Class' by LEfSe. [required]"
    )

    prepare_lefse_parser.add_argument(
        "-s",
        "--subclass-col",
        metavar="TEXT",
        help="Metadata column used as 'Subclass' by LEfSe."
    )

    prepare_lefse_parser.add_argument(
        "-u",
        "--subject-col",
        metavar="TEXT",
        help="Metadata column used as 'Subject' by LEfSe."
    )

    prepare_lefse_parser.add_argument(
        "-w",
        "--where",
        metavar="TEXT",
        help="SQLite 'WHERE' clause specifying sample metadata criteria."
    )

    prepare_lefse_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    args = parser.parse_args()
    command = args.command
    delattr(args, "command")
    commands[command](**vars(args))

if __name__ == "__main__":
    main()
