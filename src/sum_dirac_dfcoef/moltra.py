from io import TextIOWrapper
from typing import ClassVar, Dict, List

from sum_dirac_dfcoef.utils import (
    delete_dirac_input_comment_out,
    is_dirac_input_keyword,
    is_dirac_input_line_should_be_skipped,
    is_dirac_input_section,
    is_dirac_input_section_two_stars,
    is_end_dirac_input_field,
    is_start_dirac_input_field,
    space_separated_parsing_upper,
)


class MoltraInfo:
    """ Class to store MOLTRA information for the sum_dirac_dfcoef module.

    Attributes:
        is_default (bool): True if the .ACTIVE section is not found
        range_str (List[str]): List of strings in the .ACTIVE section
        range_dict (Dict[str, str]): Dictionary of the .ACTIVE section. Key: symmetry type, Value: range string
    """
    is_default: bool = True
    range_str: ClassVar[List[str]] = []
    # range_str Example:
    # ['energy -20 10 2', '10..180', ...]
    range_dict: ClassVar[Dict[str, str]] = {}

    @classmethod
    def read_moltra_section(cls, dirac_output: TextIOWrapper):
        """Read the MOLTRA section settings from the output file of DIRAC

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
            None (cls.range_str and cls.is_default will be updated)
        """

        is_moltra_section = False
        is_reach_input_field = False
        is_next_line_active = False
        for line in dirac_output:
            no_comment_line = delete_dirac_input_comment_out(line)
            words = space_separated_parsing_upper(no_comment_line)
            if is_dirac_input_line_should_be_skipped(words):
                continue

            if is_start_dirac_input_field(line):
                is_reach_input_field = True
                continue

            if is_end_dirac_input_field(line):
                break  # end of input field

            if is_reach_input_field:
                if is_dirac_input_section_two_stars(words[0]):
                    if "**MOLTRA" in words[0]:
                        is_moltra_section = True
                        continue

                if is_moltra_section:
                    if ".ACTIVE" in words[0]:
                        cls.is_default = False
                        is_next_line_active = True
                        continue

                if is_next_line_active:
                    if is_dirac_input_section(words[0]) or is_dirac_input_keyword(words[0]):
                        # End of the .ACTIVE section
                        break
                    cls.range_str.append(no_comment_line.strip())
