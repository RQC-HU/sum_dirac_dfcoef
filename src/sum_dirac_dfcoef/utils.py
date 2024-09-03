import pickle
import re
import sys
from pathlib import Path
from typing import Any, List


def space_separated_parsing(line: str) -> List[str]:
    return [word for word in line.rstrip("\n").split(" ") if word != ""]


def space_separated_parsing_upper(line: str) -> List[str]:
    return [word.upper() for word in line.rstrip("\n").split(" ") if word != ""]


def fast_deepcopy_pickle(obj: Any) -> Any:
    """Fast deep copy using pickle.

    Args:
        obj (Any): object to be copied

    Returns:
        Any: copied object
    """
    return pickle.loads(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))  # noqa: S301


def debug_print(message: str) -> None:
    from sum_dirac_dfcoef.args import args

    # print debug message if --debug option is used
    if args.debug:
        print(message)


def is_float(parameter: str) -> bool:
    try:
        float(parameter)
        return True
    except ValueError:
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

    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        sys.exit(f"ERROR: DIRAC output file is not found. file={input_path}")
    elif input_path.is_dir():
        sys.exit(
            "ERROR: The path you specified as the DIRAC output is a directory. Not a file.\
Please check your -i or --input option is correct."
        )
    return input_path


def should_write_positronic_results_to_file() -> bool:
    from sum_dirac_dfcoef.args import args

    if args.all_write or args.positronic_write:
        return True
    else:
        return False


def should_write_electronic_results_to_file() -> bool:
    from sum_dirac_dfcoef.args import args

    if args.all_write or not args.positronic_write:
        return True
    else:
        return False
