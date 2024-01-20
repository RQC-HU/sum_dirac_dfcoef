import re
from collections import OrderedDict
from enum import Enum, auto
from io import TextIOWrapper
from typing import Dict, List, Optional, Tuple
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.utils import (
    debug_print,
    delete_dirac_input_comment_out,
    is_dirac_input_keyword,
    is_dirac_input_line_should_be_skipped,
    is_dirac_input_section_one_star,
    is_end_dirac_input_field,
    is_start_dirac_input_field,
    space_separated_parsing,
    space_separated_parsing_upper,
)


class StageEigenvalues(Enum):
    INIT = auto()
    SEARCH_EIGENVALUES_HEADER = auto()
    SEARCH_PRINT_TYPE = auto()
    EIGENVALUES_READ = auto()
    OCCUPATION_INFO_READ = auto()
    WAIT_END = auto()


class StageEIGPRI(Enum):
    INIT = auto()
    SEARCH_SCF_SECTION = auto()
    SEARCH_SCF_DETAIL_SECTION = auto()
    SEARCH_EIGPRI = auto()
    NEXT_LINE_EIGPRI = auto()
    WAIT_END = auto()


# type definition eigenvalues.shell_num
# type shell_num = {
#     "E1g": {
#         "closed": int
#         "open": int
#         "virtual": int
#     },
#     "E1u": {
#         "closed": int
#         "open": int
#         "virtual": int
#     },
# }
class Eigenvalues:
    shell_num: ODict[str, Dict[str, int]]
    energies: ODict[str, Dict[int, float]]
    energies_used: ODict[str, Dict[int, bool]]

    def __init__(
        self,
        shell_num: Optional[ODict[str, Dict[str, int]]] = None,
        energies: Optional[ODict[str, Dict[int, float]]] = None,
        energies_used: Optional[ODict[str, Dict[int, bool]]] = None,
    ) -> None:
        self.shell_num = shell_num if shell_num is not None else OrderedDict()
        self.energies = energies if energies is not None else OrderedDict()
        self.energies_used = energies_used if energies_used is not None else OrderedDict()

    def __repr__(self) -> str:
        return f"shell_num: {self.shell_num}\nenergies: {self.energies}\nenergies_used: {self.energies_used}"

    def setdefault(self, key: str):
        self.shell_num.setdefault(key, {"closed": 0, "open": 0, "virtual": 0, "negative": 0, "positronic": 0})
        self.energies.setdefault(key, {})
        self.energies_used.setdefault(key, {})

    def get_electronic_spinor_num(self, symmetry_type: str) -> int:
        return self.shell_num[symmetry_type]["closed"] + self.shell_num[symmetry_type]["open"] + self.shell_num[symmetry_type]["virtual"]

    def get_eigenvalues(self, dirac_output: TextIOWrapper):
        def is_end_of_read(line) -> bool:
            return True if "HOMO - LUMO" in line else False

        def is_eigenvalue_type_written(words: List[str]) -> bool:
            # closed shell: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1043
            # open shell: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1053
            # virtual eigenvalues: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1064
            # negative energy eigenvalues (only atom or linear molecule case): https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1156
            # positronic eigenvalues (not atom and linear molecule case): https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1073
            if "*" == words[0] and "Closed" == words[1] and "shell," == words[2]:
                return True
            elif "*" == words[0] and "Open" == words[1] and "shell" == words[2]:
                return True
            elif "*" == words[0] and "Virtual" == words[1] and "eigenvalues," == words[2]:
                return True
            elif "*" == words[0] and "Negative" == words[1] and "energy" == words[2] and "eigenvalues," == words[3]:
                return True
            elif "*" == words[0] and "Positronic" == words[1] and "eigenvalues," == words[2]:
                return True
            return False

        def is_occupation_info_header(line: str) -> bool:
            return True if "Occupation in fermion symmetry" in line else False

        def prepare_occupation_info(words: List[str]) -> Tuple[str, Dict[str, int]]:
            current_symmetry_type = words[len(words) - 1]
            occ_idx = {k: 1 for k in omega[current_symmetry_type].keys()}
            return current_symmetry_type, occ_idx

        def get_current_eigenvalue_type(words: List[str]) -> str:
            # words[0] = '*', words[1] = "Closed" or "Open" or "Virtual" or "Negative" or "Positronic"
            current_eigenvalue_type = words[1].lower()
            return current_eigenvalue_type

        def get_symmetry_type_standard(words: List[str]) -> str:
            current_symmetry_type = words[3]
            return current_symmetry_type

        def get_symmetry_type_supersym(words: List[str]) -> str:
            # https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1097-1105
            # FORMAT '(/A,I4,4A,I2,...)'
            # DATA "* Block",ISUB,' in ',FREP(IFSYM),":  ",...
            # ISUB might be **** if ISUB > 9999 or ISUB < -999 because of the format
            # Therefore, find 'in' word list and get FREP(IFSYM) from the word list
            # FREP(IFSYM) is a symmetry type
            idx = words.index("in")
            current_symmetry_type = words[idx + 1][: len(words[idx + 1]) - 1]
            return current_symmetry_type

        def get_omega_str(words: List[str]) -> str:
            if "Omega" in line:
                # * Block   3 in E1u:  Omega =  5/2
                # => 5/2
                omega_str = words[len(words) - 1].replace("=", "").strip()
                return omega_str
            else:
                # * Block   3 in E1u:  p 3/2; -3/2
                # => p 3/2 -3/2
                colon_idx = line.index(":")
                omega_str = line[colon_idx + 1 : len(line) - 1].strip()
                if ";" in omega_str:
                    jval = omega_str.split(";")[0].strip()
                    mjval = omega_str.split(";")[1].strip()
                    omega_str = f"{jval} {mjval}"
                return omega_str

        def supersym_append_energies() -> None:
            for item in omega_list:
                if item not in occ_idx.keys():
                    msg = f"Cannot find {item} in occ_idx.keys()!, occ_idx.keys(): {occ_idx.keys()}"
                    raise ValueError(msg)
                val = omega[current_symmetry_type][item][occ_idx[item]]
                idx = len(self.energies[current_symmetry_type]) + 1
                self.energies[current_symmetry_type][idx] = val
                occ_idx[item] += 1

        def create_splitted_by_slash2_list(line: str) -> List[str]:
            # split by /2 and remove empty string
            # (e.g. 1, atomic  ) "s 1/2 d 3/2 s 1/2" => ["s 1/2", "d 3/2", "s 1/2"]
            # (e.g. 2, molecule) "1/2 1/2 1/2 3/2 1/2 5/2" => ["1/2", "1/2", "1/2", "3/2", "1/2", "5/2"]
            split_by_slash2 = list(filter(None, [item.strip("\r\n") for item in line.split("/2")]))
            split_by_slash2 = [f"{item.strip()}/2" for item in split_by_slash2]
            return split_by_slash2

        def add_mj_info_to_omega_list(line: str, omega_list: List[str]) -> List[str]:
            mj_list = create_splitted_by_slash2_list(line.replace("Mj", "").strip("\r\n"))
            if len(mj_list) != len(omega_list):
                msg = f"len(mj_list) != len(omega_list)\nmj_list: {mj_list}\nomega_list: {omega_list}\nline: {line}\n"
                raise ValueError(msg)
            omega_list = [f"{omega_list[i]} {mj_list[i]}" for i in range(len(mj_list))]
            return omega_list

        def read_eigenvalues(line: str) -> None:
            start_idx = 0
            while True:
                # e.g. -775.202926514  ( 2) => -775.202926514
                regex = r"[-]?[0-9]+\.?[0-9]+"
                match = re.search(regex, line[start_idx:])
                if match is None:
                    break
                val = float(match.group())

                # e.g. -775.202926514  ( 2) => 2
                regex = r"\([ ]*[0-9]+\)"
                match = re.search(regex, line[start_idx:])
                if match is None:
                    break
                # match.group() == ( 2) => [1 : len(match.group()) - 1] == 2
                num = int(match.group()[1 : len(match.group()) - 1])
                self.shell_num[current_symmetry_type][current_eigenvalue_type] += num
                if print_type == "standard":
                    for _ in range(0, num, 2):
                        idx = len(self.energies[current_symmetry_type]) + 1
                        self.energies[current_symmetry_type][idx] = val
                elif print_type == "supersymmetry":
                    for _ in range(0, num, 2):
                        idx = len(omega[current_symmetry_type][omega_str]) + 1
                        omega[current_symmetry_type][omega_str][idx] = val
                start_idx += match.end()

        stage = StageEigenvalues.INIT
        atomic = False
        print_type = ""  # "standard" or "supersymmetry"
        occ_idx: Dict[str, int] = {}
        omega: Dict[str, Dict[str, Dict[int, float]]] = {}
        omega_str = ""  # 1/2 or 3/2 or 5/2 or p 3/2 -3/2 ...
        omega_list = []
        current_eigenvalue_type = ""  # "closed" or "open" or "virtual"
        current_symmetry_type = ""  # "E1g" or "E1u" or "E1" ...
        eigenvalue_type_omega_replacement = {"inactive": "closed", "active": "open", "virtual": "virtual"}

        for line in dirac_output:
            words: List[str] = space_separated_parsing(line)

            if len(words) == 0:
                continue
            elif stage == StageEigenvalues.INIT:
                if "SCF - CYCLE" in line:
                    stage = StageEigenvalues.SEARCH_EIGENVALUES_HEADER
            elif stage == StageEigenvalues.SEARCH_EIGENVALUES_HEADER:
                if "Eigenvalues" == words[0]:
                    stage = StageEigenvalues.SEARCH_PRINT_TYPE
            elif stage == StageEigenvalues.SEARCH_PRINT_TYPE:
                if "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                    stage = StageEigenvalues.EIGENVALUES_READ
                    print_type = "standard"
                    current_symmetry_type = get_symmetry_type_standard(words)
                    self.setdefault(current_symmetry_type)
                elif "* Block" in line:
                    stage = StageEigenvalues.EIGENVALUES_READ
                    print_type = "supersymmetry"
                    current_symmetry_type = get_symmetry_type_supersym(words)
                    atomic = ";" in line
                    omega_str = get_omega_str(words)
                    self.setdefault(current_symmetry_type)
                    omega.setdefault(current_symmetry_type, {}).setdefault(omega_str, {})
            elif is_end_of_read(line) or stage == StageEigenvalues.WAIT_END:
                break
            elif stage == StageEigenvalues.EIGENVALUES_READ:
                if print_type == "standard" and "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                    current_symmetry_type = get_symmetry_type_standard(words)
                    self.setdefault(current_symmetry_type)
                elif print_type == "supersymmetry" and "* Block" in line:
                    current_symmetry_type = get_symmetry_type_supersym(words)
                    omega_str = get_omega_str(words)
                    self.setdefault(current_symmetry_type)
                    omega.setdefault(current_symmetry_type, {}).setdefault(omega_str, {})
                elif is_eigenvalue_type_written(words):
                    current_eigenvalue_type = get_current_eigenvalue_type(words)
                elif is_occupation_info_header(line):
                    stage = StageEigenvalues.OCCUPATION_INFO_READ
                    current_symmetry_type, occ_idx = prepare_occupation_info(words)
                else:
                    read_eigenvalues(line)
            elif stage == StageEigenvalues.OCCUPATION_INFO_READ:
                if "Occupation of" in line:
                    stage = StageEigenvalues.WAIT_END
                elif is_occupation_info_header(line):
                    current_symmetry_type, occ_idx = prepare_occupation_info(words)
                elif "orbitals" in line:
                    # * Inactive orbitals => inactive
                    occ_type = words[1].lower()
                    # inactive => closed
                    current_eigenvalue_type = eigenvalue_type_omega_replacement[occ_type]
                elif "Mj" in line:
                    omega_list = add_mj_info_to_omega_list(line, omega_list)
                    supersym_append_energies()
                elif atomic:
                    omega_list = create_splitted_by_slash2_list(line)
                else:  # molecular
                    omega_list = create_splitted_by_slash2_list(line)
                    supersym_append_energies()

        for key in self.energies.keys():
            num = len(self.energies[key])
            self.energies_used[key] = {i: False for i in range(1, num + 1)}

        debug_print(f"eigenvalues: {self}")

    def validate_eigpri_option(self, dirac_output: TextIOWrapper):
        """Validate the .EIGPRI option in the DIRAC input file,
        if is not set, it is a valid input
        because only the positive energy eigenvalues are printed as default.

        Args:
            dirac_output (TextIOWrapper): _description_
        """

        stage = StageEIGPRI.INIT
        for line in dirac_output:
            no_comment_out_line = delete_dirac_input_comment_out(line)
            words = space_separated_parsing_upper(no_comment_out_line)
            if is_dirac_input_line_should_be_skipped(words):
                continue

            if is_end_dirac_input_field(no_comment_out_line) or stage == StageEIGPRI.WAIT_END:
                break

            if stage == StageEIGPRI.INIT:
                if is_start_dirac_input_field(no_comment_out_line):
                    stage = StageEIGPRI.SEARCH_SCF_SECTION

            if stage == StageEIGPRI.SEARCH_SCF_SECTION:
                if is_dirac_input_keyword(words[0]) and ".SCF" in words[0]:
                    stage = StageEIGPRI.SEARCH_SCF_DETAIL_SECTION

            if stage == StageEIGPRI.SEARCH_SCF_DETAIL_SECTION:
                if is_dirac_input_section_one_star(words[0]) and "*SCF" in words[0]:
                    stage = StageEIGPRI.SEARCH_EIGPRI

            if stage == StageEIGPRI.SEARCH_EIGPRI:
                if is_dirac_input_keyword(words[0]) and ".EIGPRI" in words[0]:
                    stage = StageEIGPRI.NEXT_LINE_EIGPRI

            if stage == StageEIGPRI.NEXT_LINE_EIGPRI:
                # https://diracprogram.org/doc/master/manual/wave_function/scf.html#eigpri
                if len(words) == 2 and words[0].isdigit() and words[1].isdigit():
                    if int(words[0]) == 0:  # positive energy eigenvalues are not printed
                        msg = f"\nYour .EIGPRI option in your DIRAC input file is invalid!\n\
.EIGPRI\n\
{line}\n\
We cannot get the eigenvalues with your .EIGPRI option.\n\
If you want to use this output file with this program, you must use --no-scf option to skip reading eigenvalues information.\n\
But you cannot use the output using --no-scf option to dcaspt2_input_generator program.\n\
If you want to get the eigenvalues information, please refer the .EIGPRI option in the manual of DIRAC.\n\
https://diracprogram.org/doc/master/manual/wave_function/scf.html#eigpri\n\
You must enable to print out the positive eigenvalues energy.\n"
                        raise ValueError(msg)
                    stage = StageEIGPRI.WAIT_END
