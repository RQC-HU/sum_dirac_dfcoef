import re
from collections import OrderedDict
from io import TextIOWrapper
from typing import ClassVar, List
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.utils import (
    delete_comment_out,
    is_dirac_input_keyword,
    is_dirac_input_line_comment_out,
    is_dirac_input_section,
    is_dirac_input_section_two_stars,
    space_separated_parsing,
)


class MoltraInfo:
    is_default: bool = True
    range_str: ClassVar[List[str]] = []
    # range_str Example:
    # ['energy -20 10 2', '10..180', ...]
    idx_range: ClassVar[ODict[str, List[int]]] = OrderedDict()

    @classmethod
    def read_moltra_section(cls, dirac_output: TextIOWrapper):
        """Read the MOLTRA section settings from the output file of DIRAC

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
            None (cls.range_str and cls.is_default will be updated)
        """

        is_moltra_section = False
        re_active_keyword = r" *\.ACTIVE"
        is_next_line_active = False
        for line in dirac_output:
            words = [word.upper() for word in space_separated_parsing(line)]
            if len(words) == 0 or is_dirac_input_line_comment_out(words[0]):
                continue

            if is_dirac_input_section_two_stars(words[0]):
                if "**MOLTRA" in words[0]:
                    is_moltra_section = True
                    continue

            if is_moltra_section:
                if re.match(re_active_keyword, line) is not None:
                    cls.is_default = False
                    is_next_line_active = True
                    continue

            if is_next_line_active:
                if is_dirac_input_section(words[0]) or is_dirac_input_keyword(words[0]):
                    # End of the .ACTIVE section
                    break
                no_comment_line = delete_comment_out(line)
                cls.range_str.append(no_comment_line.strip())
