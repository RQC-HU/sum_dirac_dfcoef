import argparse
import sys


class PrintVersionExitAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):  # noqa: A002
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):  # noqa: ARG002
        from sum_dirac_dfcoef.__about__ import __version__

        print(f"{__version__}")
        sys.exit()


class PrintHelpArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stdout)
        err_msg = f"{self.prog}: error: {message}\n"
        self.exit(2, err_msg)


def parse_args() -> "argparse.Namespace":
    parser = PrintHelpArgumentParser(
        description="Summarize the coefficients from DIRAC output file that *PRIVEC option is used. (c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)"
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="(required) file name of DIRAC output", dest="file")
    parser.add_argument("-o", "--output", type=str, help="Output file name. Default: sum_dirac_dfcoef.out", dest="output")
    parser.add_argument(
        "-g",
        "--for-generator",
        action="store_true",
        help="Automatically set the arguments for dcaspt2_input_generator. \
This option is useful when you want to use the result of this program as input to dcaspt2_input_generator. \
This option is equivalent to set -c/--compress and not set -p/--positronic and --no-scf options.",
        dest="for_generator",
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        help="Compress output. Display all coefficients on one line for each MO. This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.",
        dest="compress",
    )
    parser.add_argument(
        "--only-moltra",
        action="store_true",
        help="Print only MOs that is included in the range of MOLTRA. You should activate this option when you want to get compressed output (-c/--compress option)\
            but you don't want to get the output that is not included in the range of MOLTRA.",
        dest="only_moltra",
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.1, help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution", dest="threshold"
    )
    parser.add_argument(
        "-d",
        "--decimal",
        type=int,
        default=5,
        choices=range(1, 16),
        help="Set the decimal places. Default: 5 (e.g) --decimal=3 => print orbital with 3 decimal places (0.123, 2.456, ...). range: 1-15",
        dest="decimal",
    )
    parser.add_argument("-a", "--all-write", action="store_true", help="Print all MOs(Positronic and Electronic).", dest="all_write")
    parser.add_argument(
        "-p",
        "--positronic-write",
        action="store_true",
        help="Print only Positronic MOs. The output with this option cannot be used as input to dcaspt2_input_generator.",
        dest="positronic_write",
    )
    parser.add_argument("-v", "--version", action=PrintVersionExitAction, help="Print version and exit", dest="version")
    parser.add_argument(
        "--no-scf",
        action="store_true",
        help="If you don't activate .SCF keyword in your DIRAC input file, you must use this option.\
            But you cannot use the output using this option to dcaspt2_input_generator program.",
        dest="no_scf",
    )
    parser.add_argument("--debug", action="store_true", help="print debug output (Normalization constant, Sum of MO coefficient)", dest="debug")
    parser.add_argument("--no-sort", action="store_true", help="Don't sort the output by MO energy")
    # If -v or --version option is used, print version and exit
    args = parser.parse_args()

    if args.for_generator:
        args.no_scf = False
        args.compress = True
        args.positronic_write = False

    if not (args.no_scf or args.positronic_write) and args.compress:
        args.for_generator = True

    if args.only_moltra and args.for_generator:
        parser.error("--only-moltra option cannot be used with --for-generator option.\nUse either --only-moltra or --for-generator option.")

    if args.all_write and args.positronic_write:
        parser.error("-a/--all-write and -p/--positronic-write options cannot be set at the same time.")

    if args.only_moltra and not args.compress:
        print("Warning: --only-moltra option is activated but --compress option is not activated. --only-moltra option will be ignored.")

    return args


args = parse_args()
