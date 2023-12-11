import re
from io import TextIOWrapper

from sum_dirac_dfcoef.utils import space_separated_parsing


def get_electron_num_from_input(dirac_output: TextIOWrapper) -> int:
    electron_num: int = 0
    is_next_line_electron_num: bool = False
    is_next_line_print_setting: bool = False
    is_reach_input_field: bool = False
    is_scf_found: bool = False
    scf_detail_section: bool = False
    regex_number = r"[0-9]+"
    # *section name or **section name
    regex_section = r"^\*{1,2}[0-9A-Z]+"
    for line in dirac_output:
        words = space_separated_parsing(line)
        words = [word.upper() for word in words]
        # find "Contents of the input file"
        if "Contents of the input file" in line:
            is_reach_input_field = True
            continue

        # find "Contents of the molecule file"
        if "Contents of the molecule file" in line:
            break  # end of input field

        if is_reach_input_field:
            if len(words) == 0:
                continue
            if ".SCF" in words[0]:
                is_scf_found = True

            if re.match(regex_section, words[0]) is not None:
                if "*SCF" in words[0]:
                    scf_detail_section = True
                else:
                    scf_detail_section = False

            if scf_detail_section:
                if is_next_line_electron_num:
                    number = re.search(regex_number, words[0]).group()
                    electron_num += int(number)
                    is_next_line_electron_num = False

                if is_next_line_print_setting:
                    # https://gitlab.com/kohei-noda/dirac/-/blob/79e6b9e27cf8018999ddca2aa72247ccfb9d2a2d/src/dirac/dirrdn.F#L2865-2876
                    # ipreig = 0 (no eigenvalue printout) and 2 (only positronic eigenvalues written out)
                    # are not supported because we cannot get electron number from them
                    number = re.search(regex_number, words[0]).group()
                    ipreig = int(number)
                    if ipreig in (0, 2):
                        msg = ".PRINT setting in *SCF section with value 0 or 2 is not supported.\n\
0 means no eigenvalue printout and 2 means only positronic eigenvalues written out.\n\
We cannot get electron number from them.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
    Please check your DIRAC input file and try again.\n"
                        raise ValueError(msg)
                    is_next_line_print_setting = False

                if ".PRINT" in line.upper():
                    is_next_line_print_setting = True

                # .CLOSED SHELL
                if ".CLOSED" == words[0] and "SHELL" in words[1]:
                    is_next_line_electron_num = True

                # .OPEN SHELL
                if ".OPEN" == words[0] and "SHELL" in words[1]:
                    is_next_line_electron_num = True
    if not is_scf_found:
        msg = "Cannot find SCF calculation settings from your DIRAC input file wrtte in your output file.\n\
we cannot get information about the electron number and orbital energy without SCF calculation.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
Please check your DIRAC input file and try again.\n"
        raise ValueError(msg)
    return electron_num


def get_electron_num_from_scf_field(dirac_output: TextIOWrapper) -> int:
    # https://gitlab.com/kohei-noda/dirac/-/blob/79e6b9e27cf8018999ddca2aa72247ccfb9d2a2d/src/dirac/dirrdn.F#L2127
    # find "i.e. no. of electrons ="
    is_wave_function_module_reached: bool = False
    for line in dirac_output:
        words = space_separated_parsing(line)
        if "Wave function module" in line:
            is_wave_function_module_reached = True
            continue

        if is_wave_function_module_reached:
            if "i.e. no. of electrons" in line:
                # ["i.e.", "no.", "of", "electrons", "=", number]
                return int(words[5])
    msg = "Cannot find electron number from your DIRAC output file.\n\
we cannot get information about the electron number and orbital energy without SCF calculation.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
Please check your DIRAC input file and try again.\n"
    raise ValueError(msg)
