def space_separated_parsing(line: str) -> "list[str]":
    import re

    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def debug_print(str: str):
    from .args import args

    # print debug message if --debug option is used
    if args.debug:
        print(str)
