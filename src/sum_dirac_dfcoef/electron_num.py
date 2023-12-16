import re
from io import TextIOWrapper

from sum_dirac_dfcoef.utils import is_dirac_input_section, is_end_dirac_input_field, is_start_dirac_input_field, space_separated_parsing


def get_electron_num_from_input(dirac_output: TextIOWrapper) -> int:
    """If users calculate SCF with open shell they must explicitly write the OPEN SHELL and CLOSED SHELL keywords
    in the input file. Therefore, we can get the electron number from the input file.

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
            int: The number of electrons in the system
    """

    def get_an_natural_number(word: str) -> int:
        # OPEN SHELL and CLOSED SHELL electron number settings
        # must be written as a natural number (negative number is not allowed)
        regex_number = r"[-]?[0-9]+"
        number = int(re.search(regex_number, word).group())
        if number < 0:
            msg = "The number of electrons in OPEN SHELL and CLOSED SHELL must be a natural number.\n\
But we found a negative number in your DIRAC input file.\n\
Please check your DIRAC input file and try again.\n"
            raise ValueError(msg)
        return number

    electron_num: int = 0
    is_closed_shell_section: bool = False
    is_openshell_section: bool = False
    num_of_open_shell: int = 0
    is_next_line_print_setting: bool = False
    is_reach_input_field: bool = False
    is_scf_found: bool = False
    scf_detail_section: bool = False
    # *section name or **section name
    regex_scf_keyword = r" *\.SCF"
    regex_comment_out = r" *[!#]"
    for line in dirac_output:
        words = space_separated_parsing(line)
        words = [word.upper() for word in words]

        if is_start_dirac_input_field(line):
            is_reach_input_field = True
            continue

        if is_end_dirac_input_field(line):
            break  # end of input field

        if is_reach_input_field:
            if len(words) == 0 or re.match(regex_comment_out, words[0]) is not None:
                # Comment out line or empty line
                continue

            if re.match(regex_scf_keyword, words[0]) is not None:
                is_scf_found = True

            if is_dirac_input_section(words[0]):
                if "*SCF" in words[0]:
                    scf_detail_section = True
                else:
                    scf_detail_section = False

            if scf_detail_section:
                if is_openshell_section:
                    if num_of_open_shell == 0:
                        num_of_open_shell = get_an_natural_number(words[0])
                    else:
                        num_of_open_shell -= 1
                        # open shell format
                        # https://diracprogram.org/doc/master/manual/wave_function/scf.html#open-shell
                        # .OPEN SHELL
                        # num_of_open_shell
                        # num_of_elec/irrep1_num_spinor irrep2_num_spinor ...
                        # We want to get only num_of_elec
                        electron_num += get_an_natural_number(words[0])
                        if num_of_open_shell == 0:
                            is_openshell_section = False

                if is_closed_shell_section:
                    # closed shell format
                    # https://diracprogram.org/doc/master/manual/wave_function/scf.html#closed-shell
                    # .CLOSED SHELL
                    # irrep1_num_spinor irrep2_num_spinor ...
                    # "!" or "#" means comment out(see https://gitlab.com/dirac/dirac/-/blob/ea717cdb294035d8af3ebe2b1e00cf94f1c1a6b7/src/input/parse_input.F90#L53-54)
                    # So we need to read electron numbers until the line contains "!" or "#" or the last element of words.
                    regex_comment_out = r"[!#]"
                    for word in words:
                        comment_out_str = re.search(regex_comment_out, word)
                        if comment_out_str is not None:
                            comment_idx = comment_out_str.start()
                            w = word[:comment_idx]
                            if len(w) > 0:
                                electron_num += get_an_natural_number(w)
                            break  # end of closed shell section because we found comment out
                        else:
                            electron_num += get_an_natural_number(word)
                    is_closed_shell_section = False

                if is_next_line_print_setting:
                    # https://gitlab.com/kohei-noda/dirac/-/blob/79e6b9e27cf8018999ddca2aa72247ccfb9d2a2d/src/dirac/dirrdn.F#L2865-2876
                    # ipreig = 0 (no eigenvalue printout) and 2 (only positronic eigenvalues written out)
                    # are not supported because we cannot get electron number from them
                    number = get_an_natural_number(words[0])
                    ipreig = int(number)
                    if ipreig in (0, 2):
                        msg = ".PRINT setting in *SCF section with value 0 or 2 is not supported.\n\
0 means no eigenvalue printout and 2 means only positronic eigenvalues written out.\n\
Therefore we cannot get the information (e.g. orbital energies) from DIRAC output.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
Please check your DIRAC input file and try again.\n"
                        raise ValueError(msg)
                    is_next_line_print_setting = False

                if ".PRINT" in line.upper():
                    is_next_line_print_setting = True

                # .CLOSED SHELL
                if ".CLOSED" == words[0] and "SHELL" in words[1]:
                    is_closed_shell_section = True

                # .OPEN SHELL
                if ".OPEN" == words[0] and "SHELL" in words[1]:
                    is_openshell_section = True
    if not is_scf_found:
        msg = "\nCannot find SCF calculation settings from your DIRAC input file written in your output file.\n\
we cannot get information about the electron number and orbital energy without SCF calculation.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
Please check your DIRAC input file and try again.\n\
If you want to use this program by using the output file without SCF calculation, please use --no-scf option.\n\
But you cannot use the output using --no-scf option to dcaspt2_input_generator program.\n"
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
    msg = "\nCannot find electron number from your DIRAC output file.\n\
we cannot get information about the electron number and orbital energy without SCF calculation.\n\
So we cannot continue this program because we need electron number and orbital energies to summarize DIRAC output.\n\
Please check your DIRAC input file and try again.\n\
If you want to use this program by using the output file without SCF calculation, please use --no-scf option.\n\
But you cannot use the output using --no-scf option to dcaspt2_input_generator program.\n"
    raise ValueError(msg)
