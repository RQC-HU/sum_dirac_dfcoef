import argparse


def space_separated_parsing(line: str) -> "list[str]":
    import re

    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def debug_print(args: "argparse.Namespace", str: str):
    # print debug message if --debug option is used
    if args.debug:
        print(str)
