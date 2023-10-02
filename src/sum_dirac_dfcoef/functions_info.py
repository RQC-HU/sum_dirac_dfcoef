from collections import OrderedDict
import copy
from io import TextIOWrapper
import re
from typing import OrderedDict as ODict

from .atoms import ao, is_different_atom, AtomInfo
from .utils import debug_print, space_separated_parsing


class Function:
    def __init__(self, component_func: str, symmetry: str, atom: str, gto_type: str, start_idx: int, num_functions: int, multiplicity: int) -> None:
        self.component_func = component_func  # "large" or "small"
        self.symmetry = symmetry  # e.g. "Ag"
        self.atom = atom  # e.g. "Cl"
        self.gto_type = gto_type  # e.g. "dxz"
        self.start_idx = start_idx  # e.g. 1
        self.num_functions = num_functions  # e.g. 3
        self.multiplicity = multiplicity  # e.g. 2

    def get_identifier(self) -> str:
        return f"{self.component_func} {self.symmetry} {self.atom} {self.gto_type}"


class FunctionsInfo(ODict[str, ODict[str, ODict[str, ODict[int, AtomInfo]]]]):
    # FunctionsInfo(OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[int, AtomInfo]]]]
    # "large": {
    #     "Ag": {
    #         "Cl": {
    #            "1": {
    #                 AtomInfo: {
    #                     mul: 2,
    #                     functions: {
    #                         "s": 3,
    #                         "p": 3,
    #                 }
    #            },
    #            "3": {
    #                 AtomInfo: {
    #                     mul: 2,
    #                     functions: {
    #                         "s": 3,
    #                         "p": 3,
    #                     }
    #                 },...
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

    def get_number_of_info(words: "list[str]", line_str: str) -> None:
        # Parse the line_str and get the number of orbitals
        # (e.g.) Number of orbitals in each symmetry:          108    65    65    24    78    29    29    13
        #
        def get_orbital_type() -> str:
            # orbital_type is between "Number of" and "in each symmetry:"
            start_orbital_type_idx = line_str.find("Number of") + len("Number of")
            end_orbital_type_idx = line_str.find("in each symmetry:")
            if start_orbital_type_idx == -1 or end_orbital_type_idx == -1:
                return ""  # This is not expected but skip this line_str

            return line_str[start_orbital_type_idx:end_orbital_type_idx].strip()

        if "Number of" not in line_str:
            # End "Number of" section
            number_of_section["end"] = True
            return  # Skip this line_str and start reading large

        orbital_type = get_orbital_type()
        if orbital_type == "":
            return  # This is not expected but skip this line_str

        orbitals_str_list = words[words.index("symmetry:") + 1 :]
        try:
            orbitals = [int(i) for i in orbitals_str_list]
            number_of_info[orbital_type] = sum(orbitals)
        except (ValueError, TypeError, IndexError):
            # Probably ***** is included in the line_str
            # This is expected exception
            pass

    def get_symmetry(words: "list[str]") -> str:
        symmetry = words[1]  # e.g. "Ag"
        bra_idx = symmetry.find("(")
        if bra_idx != -1:
            symmetry = symmetry[:bra_idx]
        return symmetry

    def read_func_info(words: "list[str]", line_str: str) -> Function:
        def get_last_elem() -> "tuple[int, AtomInfo]":
            # Get the last element of the OrderedDict element with the keys
            last_elem_idx, last_elem = functions_info[component_func][symmetry][atom].popitem()
            # Need to restore the last element because it was popped
            functions_info[component_func][symmetry][atom][last_elem_idx] = last_elem
            return last_elem_idx, last_elem

        def get_start_idx() -> int:
            try:
                last_elem_idx, last_elem = get_last_elem()
                start_idx = last_elem_idx + last_elem.mul
                return start_idx
            except KeyError:
                # If the start_idx does not exist, it means that this is the first element, so that the start_idx is 1
                return 1

        def parse_plabel(plabel: str) -> "tuple[str, str, str]":
            atom = plabel[:2].strip()  # e.g. "Cm" in "Cm g400"
            subshell = plabel[3].strip()  # e.g. "g" in "Cm g400"
            gto_type = plabel[3:7].strip()  # e.g. "g400" in "Cm g400"
            return atom, subshell, gto_type

        def parse_multiplicity_label(multiplicity_label: str) -> int:
            return len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4

        """Read function information from the external variables line_str and words, which are in the scope of get_functions_info() function.

        Returns:
            Function: Function information
        """

        # ref: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
        # 3X,ILAB(1,I)," functions: ",PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)
        try:
            num_functions = int(words[0])  # ILAB(1,I) (e.g.) 3
        except ValueError as e:
            # Perhaps words[0] == "******"
            raise ValueError from e  # num_functions must be integer, so raise ValueError and exit this program
        after_functions = line_str[line_str.find("functions:") + len("functions:") :].strip()  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) "Cm g400 1+2+3+4
        plabel = after_functions[:7]  # PLABEL(I,2)(6:12) (e.g.) "Cm g400"
        atom, subshell, gto_type = parse_plabel(plabel)

        multiplicity_label = after_functions[7:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
        multiplicity = parse_multiplicity_label(multiplicity_label)  # (e.g.) 4 (1+2+3+4)
        # Set the current subshell and gto_type
        ao.current_ao.set(atom, subshell, gto_type)

        function_label = symmetry + plabel.replace(" ", "")  # symmetry + PLABEL(I,2)(6:12) (e.g.) AgCms
        if is_different_atom(ao, function_label):
            # Different atom
            ao.reset()
            ao.start_idx = get_start_idx()
            # ao was reset, so set the current subshell and gto_type again
            ao.current_ao.set(atom, subshell, gto_type)
            ao.prev_ao = copy.deepcopy(ao.current_ao)

        debug_print(f"function_label: {function_label}, ao: {ao}, start_idx: {ao.start_idx}")
        ao.function_types.add(function_label)

        return Function(component_func, symmetry, atom, gto_type, ao.start_idx, num_functions, multiplicity)

    dirac_output.seek(0)  # rewind to the beginning of the file
    start_symmetry_orbitals_section = False
    number_of_section = {"start": False, "end": False}
    number_of_info: ODict[str, int] = OrderedDict()
    component_func = ""  # "large" or "small"
    symmetry = ""
    functions_info = FunctionsInfo()
    for line_str in dirac_output:
        words: "list[str]" = space_separated_parsing(line_str)
        if len(line_str) == 0:
            continue
        elif not start_symmetry_orbitals_section:
            start_symmetry_orbitals_section = is_start_symmetry_orbitals_section(words)
        elif not number_of_section["start"]:
            number_of_section["start"] = is_start_number_of_section(words)
            get_number_of_info(words, line_str)
        elif number_of_section["start"] and not number_of_section["end"]:
            get_number_of_info(words, line_str)
        elif "component functions" in line_str:
            component_func = "large" if "Large" in line_str else ("small" if "Small" in line_str else "")
        elif "Symmetry" in line_str:
            symmetry = get_symmetry(words)
        elif "functions" in line_str:
            func = read_func_info(words, line_str)
            # Create an empty dictionary if the key does not exist
            functions_info.setdefault(component_func, OrderedDict()).setdefault(symmetry, OrderedDict()).setdefault(func.atom, OrderedDict())
            # Add function information
            if func.start_idx not in functions_info[component_func][symmetry][func.atom].keys():
                label = symmetry + func.atom
                functions_info[component_func][symmetry][func.atom][func.start_idx] = AtomInfo(func.start_idx, label, func.multiplicity)
            functions_info[component_func][symmetry][func.atom][func.start_idx].add_function(func.gto_type, func.num_functions)
        elif all(char in "* \r\n" for char in line_str) and len(re.findall("[*]", line_str)) > 0:
            # all characters in line_str are * or space or line break
            break  # Stop reading symmetry orbitals

    if not start_symmetry_orbitals_section:
        raise Exception(
            "ERROR: The \"Symmetry Orbitals\" section, which is one of the essential information sections for this program, \
is not in the DIRAC output file.\n\
Please check your DIRAC output file.\n\
Perhaps you explicitly set the .PRINT option to a negative number in one of the sections?"
        )

    return functions_info
