import re
from io import TextIOWrapper
from typing import Dict

from .utils import debug_print, space_separated_parsing


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
    scf_cycle = False
    find_eigenvalues = False
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

        if scf_cycle and not find_eigenvalues and "Eigenvalues" == words[0]:
            find_eigenvalues = True
            continue

        if find_eigenvalues and print_type == "":
            if "---" in line or len(words) == 0:
                continue
            # * Fermion symmetry
            elif "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                print_type = "standard"
                current_symmetry_type = words[3]
                eigenvalues[current_symmetry_type] = {'closed': 0, 'open': 0, 'virtual': 0}
                # read
                continue
            else:
                print_type = "supersymmetry"
                # read '* Occupation in fermion symmetry ',FREP(IFSYM)
                continue

        if print_type == "standard":
            if "*" == words[0] and "Fermion" in words[1] and "symmetry" in words[2]:
                current_symmetry_type = words[3]
                eigenvalues[current_symmetry_type] = {'closed': 0, 'open': 0, 'virtual': 0}
            elif "*" == words[0] and "Closed" in words[1] and "shell" in words[2]:
                current_eigenvalue_type = "closed"
            elif "open" in words[1] and "shell" in words[2]:
                current_eigenvalue_type = "open"
            elif "Virtual" in words[1] and "eigenvalues" in words[2]:
                current_eigenvalue_type = "virtual"
            elif current_eigenvalue_type == "":
                # Skip
                continue
            elif "HOMO - LUMO" in line:
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

        # elif print_type == "supersymmetry":
    print(f"eigenvalues: {eigenvalues}")
