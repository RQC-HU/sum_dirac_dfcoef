from enum import Enum, auto
from io import TextIOWrapper

from sum_dirac_dfcoef.utils import (
    delete_dirac_input_comment_out,
    is_dirac_input_line_should_be_skipped,
    is_dirac_input_section_two_stars,
    is_end_dirac_input_field,
    is_start_dirac_input_field,
    space_separated_parsing_upper,
)


class SchemeReaderStage(Enum):
    INIT = auto()
    SEARCH_MOLTRA_SECTION = auto()
    SEARCH_SCHEME_KEYWORD = auto()
    SCHEME_READ = auto()
    WAIT_END = auto()


class Scheme:
    """Class to store **MOLTRA > .SCHEME information for the sum_dirac_dfcoef module.

    Attributes:
        value (int): 0 if .SCHEME was not set, otherwise equals to the value of .SCHEME option.
    """

    def __init__(self) -> None:
        self.value: int = 0

    def get_scheme_num_from_input(self, dirac_output: TextIOWrapper) -> None:
        """If user explicitly set **MOLTRA > .SCHEME option (https://diracprogram.org/doc/release-23/manual/moltra.html#scheme),
        read the option value and store it to self.value.
        If this option is not specified, self.value is kept at 0 (means scheme was not set).

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
            None (self.value will be updated)
        """
        stage = SchemeReaderStage.INIT
        for line in dirac_output:
            no_comment_out_line = delete_dirac_input_comment_out(line)
            words = space_separated_parsing_upper(no_comment_out_line)

            if is_dirac_input_line_should_be_skipped(words):
                continue

            if is_end_dirac_input_field(no_comment_out_line) or stage == SchemeReaderStage.WAIT_END:
                break  # end of input field

            if stage == SchemeReaderStage.INIT:
                if is_start_dirac_input_field(no_comment_out_line):
                    stage = SchemeReaderStage.SEARCH_MOLTRA_SECTION

            elif stage == SchemeReaderStage.SEARCH_MOLTRA_SECTION:
                if is_dirac_input_section_two_stars(words[0]):
                    if "**MOLTRA" in words[0]:
                        stage = SchemeReaderStage.SEARCH_SCHEME_KEYWORD

            elif stage == SchemeReaderStage.SEARCH_SCHEME_KEYWORD:
                if ".SCHEME" in words[0]:
                    stage = SchemeReaderStage.SCHEME_READ

            elif stage == SchemeReaderStage.SCHEME_READ:
                self.value = int(words[0])
                if self.value <= 0:
                    msg = f"Invalid **MOLTRA > .SCHEME value. This option must be larger than 1, but actual .SCHEME value is {self.value}."
                    raise ValueError(msg)
                stage = SchemeReaderStage.WAIT_END
