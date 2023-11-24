import re
from io import TextIOWrapper
from typing import Dict

from .utils import space_separated_parsing


# type definition eigenvalues
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
class Eigenvalues(Dict[str, Dict[str, int]]):
    pass


def get_eigenvalues(dirac_output: TextIOWrapper):
    def is_end_of_read(line) -> bool:
        if "Occupation" in line or "HOMO - LUMO" in line:
            return True
        return False

    def is_eigenvalue_type_written(words: "list[str]") -> bool:
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

    def get_current_eigenvalue_type(words: "list[str]") -> str:
        # words[0] = '*', words[1] = "Closed" or "Open" or "Virtual"
        current_eigenvalue_type = words[1].lower()
        return current_eigenvalue_type

    def get_symmetry_type_standard(words: "list[str]") -> str:
        current_symmetry_type = words[3]
        return current_symmetry_type

    def get_symmetry_type_supersym(words: "list[str]") -> str:
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
    print_type = ""  # 'standard' or 'supersymmetry'
    eigenvalues = Eigenvalues()
    current_eigenvalue_type = ""  # 'closed' or 'open' or 'virtual'
    current_symmetry_type = ""  # 'E1g' or 'E1u' or "E1" ...

    for line in dirac_output:
        words: "list[str]" = space_separated_parsing(line)

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
                eigenvalues.setdefault(current_symmetry_type, {'closed': 0, 'open': 0, 'virtual': 0})
            elif "* Block" in line:
                print_type = "supersymmetry"
                current_symmetry_type = get_symmetry_type_supersym(words)
                eigenvalues.setdefault(current_symmetry_type, {'closed': 0, 'open': 0, 'virtual': 0})
            continue

        if print_type == "standard" and "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
            current_symmetry_type = get_symmetry_type_standard(words)
            eigenvalues.setdefault(current_symmetry_type, {'closed': 0, 'open': 0, 'virtual': 0})
        elif print_type == "supersymmetry" and "* Block" in line:
            current_symmetry_type = get_symmetry_type_supersym(words)
            eigenvalues.setdefault(current_symmetry_type, {'closed': 0, 'open': 0, 'virtual': 0})
        elif is_eigenvalue_type_written(words):
            current_eigenvalue_type = get_current_eigenvalue_type(words)
        elif is_end_of_read(line):
            break
        else:
            start_idx = 0
            while True:
                # e.g. -775.202926514  ( 2) => 2
                regex = r"\([ ]*[0-9]+\)"
                match = re.search(regex, line[start_idx:])
                if match is None:
                    break
                num = int(match.group()[1 : len(match.group()) - 1])
                eigenvalues[current_symmetry_type][current_eigenvalue_type] += num
                start_idx += match.end()

    print(f"eigenvalues: {eigenvalues}")
