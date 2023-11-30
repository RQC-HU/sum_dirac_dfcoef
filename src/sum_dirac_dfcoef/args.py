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

    def __call__(self):
        from sum_dirac_dfcoef.__about__ import __version__

        print(f"{__version__}")
        sys.exit()


def parse_args() -> "argparse.Namespace":
    parser = argparse.ArgumentParser(
        description="Summarize the coefficients from DIRAC output file that *PRIVEC option is used. (c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)"
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="(required) file name of DIRAC output", dest="file")
    parser.add_argument("-o", "--output", type=str, help="Output file name. Default: sum_dirac_dfcoef.out", dest="output")
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        help="Compress output. Display all coefficients on one line for each MO. This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.",
        dest="compress",
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
    parser.add_argument("-p", "--positronic-write", action="store_true", help="Print only Positronic MOs.", dest="positronic_write")
    parser.add_argument("-v", "--version", action=PrintVersionExitAction, help="Print version and exit", dest="version")
    parser.add_argument("--debug", action="store_true", help="print debug output (Normalization constant, Sum of MO coefficient)", dest="debug")
    parser.add_argument("--no-sort", action="store_true", help="Don't sort the output by MO energy")
    # If -v or --version option is used, print version and exit
    return parser.parse_args()


args = parse_args()
