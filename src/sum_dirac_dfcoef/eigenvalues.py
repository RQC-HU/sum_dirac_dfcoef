import re
from collections import OrderedDict
from io import TextIOWrapper
from typing import ClassVar, Dict, List
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.utils import debug_print, space_separated_parsing


# type definition eigenvalues.shell_num
# type eigenvalues = {
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
    shell_num: ClassVar[ODict[str, Dict[str, int]]] = OrderedDict()
    energies: ClassVar[ODict[str, List[float]]] = OrderedDict()

    def setdefault(self, key: str):
        self.shell_num.setdefault(key, {"closed": 0, "open": 0, "virtual": 0})
        self.energies.setdefault(key, [])

    def get_eigenvalues(self, dirac_output: TextIOWrapper):
        def is_end_of_read(line) -> bool:
            if "Occupation" in line or "HOMO - LUMO" in line:
                return True
            return False

        def is_eigenvalue_type_written(words: List[str]) -> bool:
            # closed shell: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1043
            # open shell: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1053
            # virtual eigenvalues: https://gitlab.com/dirac/dirac/-/blob/364663fd2bcc419e41ad01703fd782889435b576/src/dirac/dirout.F#L1064
            if "*" == words[0] and "Closed" == words[1] and "shell," == words[2]:
                return True
            elif "*" == words[0] and "Open" == words[1] and "shell" == words[2]:
                return True
            elif "*" == words[0] and "Virtual" == words[1] and "eigenvalues," == words[2]:
                return True
            return False

        def get_current_eigenvalue_type(words: List[str]) -> str:
            # words[0] = '*', words[1] = "Closed" or "Open" or "Virtual"
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

        scf_cycle = False
        eigenvalues_header = False
        print_type = ""  # "standard" or "supersymmetry"
        current_eigenvalue_type = ""  # "closed" or "open" or "virtual"
        current_symmetry_type = ""  # "E1g" or "E1u" or "E1" ...

        for line in dirac_output:
            words: List[str] = space_separated_parsing(line)

            if len(words) == 0:
                continue

            if "SCF - CYCLE" in line:
                scf_cycle = True
                continue

            if scf_cycle and not eigenvalues_header:
                if "Eigenvalues" == words[0]:
                    eigenvalues_header = True
                continue

            if print_type == "":  # search print type (standard or supersymmetry)
                if "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                    print_type = "standard"
                    current_symmetry_type = get_symmetry_type_standard(words)
                    self.setdefault(current_symmetry_type)
                elif "* Block" in line:
                    print_type = "supersymmetry"
                    current_symmetry_type = get_symmetry_type_supersym(words)
                    self.setdefault(current_symmetry_type)
                continue

            if print_type == "standard" and "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                current_symmetry_type = get_symmetry_type_standard(words)
                self.setdefault(current_symmetry_type)
            elif print_type == "supersymmetry" and "* Block" in line:
                current_symmetry_type = get_symmetry_type_supersym(words)
                self.setdefault(current_symmetry_type)
            elif is_eigenvalue_type_written(words):
                current_eigenvalue_type = get_current_eigenvalue_type(words)
            elif is_end_of_read(line):
                break
            else:
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
                    # [1 : len(match.group()) - 1] => ( 2) => 2
                    num = int(match.group()[1 : len(match.group()) - 1])
                    self.shell_num[current_symmetry_type][current_eigenvalue_type] += num
                    for _ in range(0, num, 2):
                        self.energies[current_symmetry_type].append(val)
                    start_idx += match.end()

        for key in self.energies.keys():
            self.energies[key].sort()
        debug_print(f"eigenvalues: {self}")
