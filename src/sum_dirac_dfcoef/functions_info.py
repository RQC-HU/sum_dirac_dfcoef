from collections import OrderedDict
import copy
import json
from io import TextIOWrapper
import re


from .utils import space_separated_parsing
from .atoms import AtomicOrbitals, is_different_atom


class Function:
    def __init__(self, component_func: str, symmetry: str, atom: str, gto_type: str, start_idx: int, num_functions: int, multiplicity: int) -> None:
        self.component_func = component_func  # "large orbitals" or "small orbitals"
        self.symmetry = symmetry  # e.g. "Ag"
        self.atom = atom  # e.g. "Cl"
        self.gto_type = gto_type  # e.g. "dxz"
        self.start_idx = start_idx  # e.g. 1
        self.num_functions = num_functions  # e.g. 3
        self.multiplicity = multiplicity  # e.g. 2

    def get_identifier(self) -> str:
        return f"{self.component_func} {self.symmetry} {self.atom} {self.gto_type}"


class AtomInfo():
    mul: int
    functions: OrderedDict[str, int]

    def __init__(self, mul: int) -> None:
        self.mul = mul
        self.functions = OrderedDict()

    def add_function(self, function: Function) -> None:
        self.functions[function.gto_type] = function.num_functions

class FunctionsInfo(OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[int, AtomInfo]]]]):
    # FunctionsInfo(OrderedDict[str, OrderedDict[str, OrderedDict[str, OrderedDict[int, OrderedDict[str, OrderedDict[str, int]]]]]]
    # "large orbitals": {
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

class AtomInfoEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, AtomInfo):
            return {
                'mul': obj.mul,
                'functions': dict(obj.functions)
            }
        return super().default(obj)

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

    def read_func_info(words: "list[str]", line_str: str) -> Function:
        def get_start_idx() -> int:
            try:
                # Get the last element of the OrderedDict element with the keys of args
                start_idx = list(functions_info[component_func][symmetry][atom].keys())[-1]
                last_elem = functions_info[component_func][symmetry][atom][start_idx]
                start_idx = start_idx + last_elem.mul
                return start_idx
            except KeyError:
                # If the start_idx does not exist, it means that this is the first element, so that the start_idx is 1
                return 1

        def read_plabel(plabel: str) -> tuple[str, str, str]:
            atom = plabel[:2].strip()  # e.g. "Cm" in "Cm g400"
            subshell = plabel[3].strip()  # e.g. "g" in "Cm g400"
            gto_type = plabel[3:7].strip()  # e.g. "g400" in "Cm g400"
            return atom, subshell, gto_type

        """Read function information from the external variables line_str and words, which are in the scope of get_functions_info() function.

        Returns:
            Function: Function information
        """

        nonlocal ao

        # ref: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
        try:
            num_functions = int(words[0])  # ILAB(1,I) (e.g.) 3
        except ValueError as e:
            # Perhaps words[0] == "******"
            raise ValueError from e  # num_functions must be integer, so raise ValueError and exit this program
        after_functions = line_str[line_str.find("functions:") + len("functions:") :].strip()  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) "Cm g400 1+2+3+4
        plabel = after_functions[:7]  # PLABEL(I,2)(6:12) (e.g.) "Cm g400"
        multiplicity_label = after_functions[7:]  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4

        atom, subshell, gto_type = read_plabel(plabel)

        # Set the current subshell and gto_type
        ao.current_ao.set(atom, subshell, gto_type)
        # if len(ao.function_types) == 0:
        #     # First function
        #     print("first function")
        #     ao.prev_ao = copy.deepcopy(ao.current_ao)

        function_label = symmetry + plabel.replace(" ", "")  # symmetry + PLABEL(I,2)(6:12) (e.g.) AgCms
        if is_different_atom(ao, function_label):
            # Different atom
            ao.reset()
            ao.start_idx = get_start_idx()
            ao.current_ao.set(atom, subshell, gto_type)
            ao.prev_ao = copy.deepcopy(ao.current_ao)
            

        print(f"function_label: {function_label}, ao: {ao}, start_idx: {ao.start_idx}")
        ao.function_types.add(function_label)

        multiplicity_label = after_functions[7:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
        multiplicity = len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4

        # start_idx = get_start_idx()
        return Function(component_func, symmetry, atom, gto_type, ao.start_idx, num_functions, multiplicity)

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
                symmetry = symmetry[:bra_idx]  # e.g. "Ag" in "Ag(1)"
        elif "functions" in line_str:
            func = read_func_info(words, line_str)
            # Create an empty dictionary if the key does not exist
            functions_info[component_func].setdefault(symmetry, OrderedDict()).setdefault(func.atom, OrderedDict())
            # Add function information
            if func.start_idx not in functions_info[component_func][symmetry][func.atom].keys():
                functions_info[component_func][symmetry][func.atom][func.start_idx] = AtomInfo(func.multiplicity)
            # functions_info[component_func][symmetry][func.atom][func.start_idx].add_function(func.gto_type, func.num_functions)
            functions_info[component_func][symmetry][func.atom][func.start_idx].add_function(func)
        # all characters in line_str are * or space or line break
        elif all(char in "* \r\n" for char in line_str) and len(re.findall("[*]", line_str)) > 0:
            break  # Stop reading symmetry orbitals

    if not start_symmetry_orbitals_section:
        raise Exception(
            "ERROR: The \"Symmetry Orbitals\" section, which is one of the essential information sections for this program, \
is not in the DIRAC output file.\n\
Please check your DIRAC output file.\n\
Perhaps you explicitly set the .PRINT option to a negative number in one of the sections?"
        )
    json.dump(functions_info, open("functions_info.json", "w"), indent=4, cls=AtomInfoEncoder)

    return functions_info
