from collections import OrderedDict
import json
from io import TextIOWrapper
from typing import Annotated
import re

from annotated_types import MaxLen
from pydantic import BaseModel


from .utils import space_separated_parsing


class AtomicOrbitals(BaseModel, validate_assignment=True):
    prev_subshell: Annotated[str, MaxLen(max_length=1)] = "s"
    current_subshell: Annotated[str, MaxLen(max_length=1)] = "s"
    function_types: "set[str]" = set()

    def reset(self):
        self.prev_subshell = "s"
        self.current_subshell = "s"
        self.function_types.clear()


class FunctionsInfo(OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[int, OrderedDict[str, int]]]]]]):
    # functions_info: OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[int, OrderedDict[str, int]]]]]]
    # "large orbitals": {
    #     "Ag": {
    #         "Cl": {
    #             "s": {
    #                 1: {"start_idx": 1, "functions": 3, "mul": 2},
    #                 2: {"start_idx": 4, "functions": 3, "mul": 2},
    #                 3: {"start_idx": 7, "functions": 3, "mul": 2},
    #             }
    #         }
    #     }
    # }
    pass


def get_functions_info(dirac_output: TextIOWrapper) -> FunctionsInfo:
    def is_start_symmetry_orbitals_section(words: "list[str]") -> bool:
        # ref: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/dirac/dirtra.F#L3654
        if len(words) == 2 and words[0] == "Symmetry" and words[1] == "Orbitals":
            return True
        return False

    def is_start_number_of_section(words: "list[str]") -> bool:
        if number_of_section["start"]:
            return False
        elif len(words) >= 6 and words[0] == "Number" and words[1] == "of" and words[2] == "orbitals" and words[3] == "in" and words[4] == "each" and words[5] == "symmetry:":
            return True
        return False

    def get_number_of_info(orbital_type: str, orbitals_str_list: "list[str]") -> None:
        try:
            orbitals = [int(i) for i in orbitals_str_list]
            number_of_info[orbital_type] = sum(orbitals)
        except (ValueError, TypeError, IndexError):
            # Probably ***** is included in the line_str
            # This is expected exception
            pass

    def is_reverse_subshell() -> bool:
        order_of_subshell = "spdfghiklmnoqrtuvwxyz"
        if order_of_subshell.index(ao.prev_subshell) > order_of_subshell.index(ao.current_subshell):
            return True
        return False

    dirac_output.seek(0)  # rewind to the beginning of the file
    start_symmetry_orbitals_section = False
    number_of_section = {"start": False, "end": False}
    number_of_info = {"orbitals": 0, "large orbitals": 0, "small orbitals": 0}
    component_func = ""  # "large orbitals" or "small orbitals"
    # functions_info = {"large orbitals": {"Ag": {"labels: {"C  s": {"functions": 3, "mul": 2}, "C  p": {"functions": 3, "mul": 2}, ...}, "check": False}, "B1g": {...}, ...}, "small orbitals": {...}}
    functions_info = FunctionsInfo()
    functions_info["large orbitals"] = OrderedDict()
    functions_info["small orbitals"] = OrderedDict()
    symmetry = ""
    ao = AtomicOrbitals()
    for line_str in dirac_output:
        words: "list[str]" = space_separated_parsing(line_str)
        if len(line_str) == 0:
            continue
        elif not start_symmetry_orbitals_section:
            start_symmetry_orbitals_section = is_start_symmetry_orbitals_section(words)
        elif not number_of_section["start"]:
            if is_start_number_of_section(words):
                number_of_section["start"] = True
                get_number_of_info("orbitals", words[6:])
        elif number_of_section["start"] and not number_of_section["end"]:
            if "Number of" not in line_str:
                # End of number of section
                number_of_section["end"] = True
                continue  # Skip this line_str and start reading large orbitals
            elif len(words) <= 7:
                continue
            orbital_type = "large orbitals" if "large" in line_str else ("small orbitals" if "small" in line_str else "")
            if orbital_type == "":
                continue  # This is not expected but skip this line_str
            get_number_of_info(orbital_type, words[7:])
        elif "component functions" in line_str:
            component_func = "large orbitals" if "Large" in line_str else ("small orbitals" if "Small" in line_str else "")
        elif "Symmetry" in line_str:
            symmetry = words[1]  # e.g. "Ag"
            bra_idx = symmetry.find("(")
            if bra_idx != -1:
                symmetry = symmetry[:bra_idx]
        elif "functions" in line_str:
            # ref: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
            try:
                num_functions = int(words[0])  # ILAB(1,I) (e.g.) 3
            except ValueError as e:
                # Perhaps words[0] == "******"
                raise ValueError from e  # num_functions must be integer
            after_functions = line_str[line_str.find("functions:") + len("functions:") :].strip()  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)
            ao.current_subshell = after_functions[3]  # e.g. "g" in "Cm g400"
            gto_type = after_functions[3:7].strip()  # e.g. "g400" in "Cm g400"
            atom = after_functions[0:2].strip()  # e.g. "Cm" in "Cm g400"
            function_label = symmetry + after_functions[:7].strip().replace(" ", "")  # symmetry + PLABEL(I,2)(6:12) (e.g.) AgCms
            if function_label in ao.function_types or is_reverse_subshell():
                # Different atom
                ao.reset()
            multiplicity_label = after_functions[7:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
            multiplicity = len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4
            ao.function_types.add(function_label)
            if symmetry not in functions_info[component_func]:
                functions_info[component_func][symmetry] = OrderedDict()
            if atom not in functions_info[component_func][symmetry]:
                functions_info[component_func][symmetry][atom] = OrderedDict()
            if gto_type not in functions_info[component_func][symmetry][atom]:
                functions_info[component_func][symmetry][atom][gto_type] = OrderedDict()
            cnt = len(functions_info[component_func][symmetry][atom][gto_type]) + 1
            start_idx = (
                1
                if len(functions_info[component_func][symmetry][atom][gto_type]) == 0
                else functions_info[component_func][symmetry][atom][gto_type][cnt - 1]["start_idx"] + functions_info[component_func][symmetry][atom][gto_type][cnt - 1]["mul"]
            )
            functions_info[component_func][symmetry][atom][gto_type][cnt] = OrderedDict({"start_idx": start_idx, "functions": num_functions, "mul": multiplicity})
            # functions_info[component_func][symmetry][function_label] = {"functions": num_functions, "mul": multiplicity}
        # all characters in line_str are * or space or line break
        elif all(char in "* \r\n" for char in line_str) and len(re.findall("[*]", line_str)) > 0:
            break  # Stop reading symmetry orbitals

    if not start_symmetry_orbitals_section:
        raise Exception(
            "ERROR: The \"Symmetry Orbitals\" section, which is one of the essential information sections for this program,\
is not in the DIRAC output file.\n\
Please check your DIRAC output file.\n\
Perhaps you explicitly set the .PRINT option to a negative number in one of the sections?"
        )
    print(f"functions_info: {functions_info}")
    json.dump(functions_info, open("functions_info.json", "w"), indent=4)

    return functions_info
