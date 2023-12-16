import re
import sys
from pathlib import Path
from typing import List


def space_separated_parsing(line: str) -> List[str]:
    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def space_separated_parsing_upper(line: str) -> List[str]:
    words = re.split(" +", line.rstrip("\n"))
    return [word.upper() for word in words if word != ""]


def debug_print(message: str):
    from sum_dirac_dfcoef.args import args

    # print debug message if --debug option is used
    if args.debug:
        print(message)


def is_float(parameter: str):
    if not parameter.isdecimal():
        try:
            float(parameter)
            return True
        except ValueError:
            return False
    else:
        return False


def is_start_dirac_input_field(line: str) -> bool:
    return True if "Contents of the input file" in line else False


def is_end_dirac_input_field(line: str) -> bool:
    return True if "Contents of the molecule file" in line else False


def is_dirac_input_keyword(word: str) -> bool:
    regex_keyword = r" *\.[0-9A-Z]+"
    return re.match(regex_keyword, word) is not None


def is_dirac_input_section(word: str) -> bool:
    regex_section = r" *\*{1,2}[0-9A-Z]+"
    return re.match(regex_section, word) is not None


def is_dirac_input_section_one_star(word: str) -> bool:
    regex_section = r" *\*[0-9A-Z]+"
    return re.match(regex_section, word) is not None


def is_dirac_input_section_two_stars(word: str) -> bool:
    regex_section = r" *\*{2}[0-9A-Z]+"
    return re.match(regex_section, word) is not None


def is_dirac_input_line_comment_out(word: str) -> bool:
    regex_comment_out = r" *[!#]"
    return re.match(regex_comment_out, word) is not None


def is_dirac_input_line_should_be_skipped(words: List[str]) -> bool:
    if len(words) == 0:
        return True
    if is_dirac_input_line_comment_out(words[0]):
        return True
    return False


def delete_dirac_input_comment_out(line: str) -> str:
    regex_comment_out = r" *[!#]"
    idx_comment_out = re.search(regex_comment_out, line)
    if idx_comment_out is None:
        return line
    return line[: idx_comment_out.start()]


def get_dirac_filepath() -> Path:
    from sum_dirac_dfcoef.args import args

    if not args.file:
        sys.exit("ERROR: DIRAC output file is not given. Please use -f option.")
    elif not Path(args.file).exists():
        sys.exit(f"ERROR: DIRAC output file is not found. file={args.file}")
    path = Path(args.file)
    return path
