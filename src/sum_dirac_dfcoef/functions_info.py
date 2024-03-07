import copy
import re
from collections import OrderedDict
from io import TextIOWrapper
from typing import List, Tuple
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.atoms import AtomInfo, FuncIndices, ao, is_different_atom
from sum_dirac_dfcoef.utils import debug_print, space_separated_parsing


class Function:
    def __init__(self, component_func: str, symmetry: str, atom: str, gto_type: str, idx_within_same_atom: int, num_functions: int, multiplicity: int) -> None:
        self.component_func = component_func  # "large" or "small"
        self.symmetry = symmetry  # e.g. "Ag"
        self.atom = atom  # e.g. "Cl"
        self.gto_type = gto_type  # e.g. "dxz"
        self.idx_within_same_atom = idx_within_same_atom  # e.g. 1
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
    #                     },
    #                 }
    #            },
    #            "3": {
    #                 AtomInfo: {
    #                     mul: 2,
    #                     functions: {
    #                         "s": 3,
    #                         "p": 3,
    #                     },
    #                 },...
    #             }
    #         }
    #     }
    # }
    pass


class SymmetryOrbitalsSummary:
    all_sym: Tuple[int]
    large_sym: Tuple[int]
    small_sym: Tuple[int]

    def __init__(self):
        self.all_sym = ()
        self.large_sym = ()
        self.small_sym = ()

    def __repr__(self):
        return f"all_sym: {self.all_sym}, large_sym: {self.large_sym}, small_sym: {self.small_sym}"


class FuncIndicesDIRAC19:
    cur_func_idx: int
    last_func_idx: int

    def __init__(self):
        self.cur_func_idx = 0
        self.last_func_idx = 0

    def __repr__(self):
        return f"cur_func_idx: {self.cur_func_idx}, last_func_idx: {self.last_func_idx}"

    def add_num_functions(self, num_functions: int) -> None:
        self.cur_func_idx += num_functions

    def set_indices(self, start: int, last: int) -> None:
        self.cur_func_idx = start
        self.last_func_idx = last

    def check_prev_sym_idx(self) -> None:
        # Before reading the next symmetry, the current symmetry index must be the same as the last symmetry index
        if self.cur_func_idx != self.last_func_idx:
            msg = f"error: The current symmetry index {self.cur_func_idx} must be the same as the last symmetry index {self.last_func_idx}"
            raise ValueError(msg)


class FuncNumDIRAC21:
    gerade: int  # Number of functions with gerade symmetry
    ungerade: int  # Number of functions with ungerade symmetry
    no_inv_sym: int  # Number of functions with no inversion symmetry

    def __init__(self):
        self.gerade = 0
        self.ungerade = 0
        self.no_inv_sym = 0

    def __repr__(self):
        return f"gerade: {self.gerade}, ungerade: {self.ungerade}, sum_all: {self.sum_all}"

    def add_num_functions(self, symmetry: str, num_functions: int) -> None:
        if symmetry[-1].lower() == "g":
            self.gerade += num_functions
        elif symmetry[-1].lower() == "u":
            self.ungerade += num_functions
        else:
            self.no_inv_sym += num_functions

    def get_current_func_idx(self, symmetry: str) -> int:
        if symmetry[-1].lower() == "g":
            return self.gerade
        elif symmetry[-1].lower() == "u":
            return self.ungerade
        else:
            return self.no_inv_sym


class FuncNumSummary:
    dirac19: FuncIndicesDIRAC19
    dirac21: FuncNumDIRAC21

    def __init__(self):
        self.dirac19 = FuncIndicesDIRAC19()
        self.dirac21 = FuncNumDIRAC21()

    def __repr__(self):
        return f"dirac19: {self.dirac19}, dirac21: {self.dirac21}"


def get_functions_info(dirac_output: TextIOWrapper) -> FunctionsInfo:
    def is_start_symmetry_orbitals_section(words: List[str]) -> bool:
        # ref: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/dirac/dirtra.F#L3654
        if len(words) == 2 and words[0] == "Symmetry" and words[1] == "Orbitals":
            return True
        return False

    def get_component_func(line_str: str) -> str:
        if "Large" in line_str:
            return "large"
        elif "Small" in line_str:
            return "small"
        else:
            msg = f"error: Unknown/Unsupported component function: {line_str}"
            raise ValueError(msg)

    def get_symmetry(words: List[str]) -> str:
        symmetry = words[1]  # e.g. "Ag"
        bra_idx = symmetry.find("(")
        if bra_idx != -1:
            symmetry = symmetry[:bra_idx]
        return symmetry

    def get_symmetry_idx(line_str: str) -> int:
        bra_idx = line_str.find("(")
        ket_idx = line_str.find(")")
        if bra_idx == -1 or ket_idx == -1:
            msg = f"error: The symmetry index is not found in the line: {line_str}"
            raise ValueError(msg)
        return int(line_str[bra_idx + 1 : ket_idx])

    def check_symmetry_idx(idx_symmetry: int, symmetry: str, summary: SymmetryOrbitalsSummary) -> None:
        if idx_symmetry < 0:
            msg = f"error: The symmetry index must be >= 0, but got {idx_symmetry} in symmetry: {symmetry}"
            raise ValueError(msg)
        elif idx_symmetry > len(summary.all_sym):
            msg = f"error: The symmetry index must be <= {len(summary.all_sym)}, but got {idx_symmetry} in symmetry: {symmetry}"
            raise ValueError(msg)

    def read_func_info(words: List[str], line_str: str, component_func: str, symmetry: str) -> Function:
        """Read function information from the external variables line_str and words
        format url: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
        actual format: 3X,ILAB(1,I)," functions: ",PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)

        (e.g.)  line_str = "18 functions:    Ar s", component_func = "large", symmetry = "Ag"
                => return Function(component_func="large", symmetry="Ag", atom="Ar", gto_type="s", idx_within_same_atom=1, num_functions=18, multiplicity=1)

                line_str = "6 functions:    H  s   1+2", component_func = "large", symmetry = "A1"
                => return Function(component_func="large", symmetry="A1", atom="H", gto_type="s", idx_within_same_atom=1, num_functions=6, multiplicity=2)

        Returns:
            Function: Function information
        """

        def get_last_elem() -> Tuple[int, AtomInfo]:
            # Get the last element of the OrderedDict element with the keys
            last_elem_atom_idx, last_elem = functions_info[component_func][symmetry][atom].popitem()
            # Need to restore the last element because it was popped
            functions_info[component_func][symmetry][atom][last_elem_atom_idx] = last_elem
            return last_elem_atom_idx, last_elem

        def get_idx_within_the_same_atom() -> int:
            try:
                last_elem_atom_idx, last_elem = get_last_elem()
                idx_within_same_atom = last_elem_atom_idx + last_elem.mul
                return idx_within_same_atom
            except KeyError:
                # If the idx_within_same_atom does not exist, it means that this is the first element, so that the idx_within_same_atom is 1
                return 1

        def parse_plabel(plabel: str) -> Tuple[str, str, str]:
            atom = plabel[4:6].strip()  # e.g. "Cm" in "    Cm g400"
            gto_type = plabel[7:].strip()  # e.g. "g400" in "    Cm g400"
            subshell = gto_type[0]  # e.g. "g" in "g400"
            return atom, subshell, gto_type

        def parse_multiplicity_label(multiplicity_label: str) -> int:
            return len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4

        try:
            num_functions = int(words[0])  # ILAB(1,I) (e.g.) 3
        except ValueError as e:
            # Perhaps words[0] == "******"
            raise ValueError from e  # num_functions must be integer, so raise ValueError and exit this program
        after_functions = line_str[line_str.find("functions:") + len("functions:") :]  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) "Cm g400 1+2+3+4
        plabel = after_functions[:11]  # "functions: ",3X,PLABEL(I,2)(6:12) (e.g.) "functions:    Cm g400" => "    Cm g400"
        atom, subshell, gto_type = parse_plabel(plabel)

        multiplicity_label = after_functions[11:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
        multiplicity = parse_multiplicity_label(multiplicity_label)  # (e.g.) 4 (1+2+3+4)
        # Set the current subshell and gto_type
        ao.current_ao.update(atom, subshell, gto_type)

        function_label = symmetry + plabel.replace(" ", "")  # symmetry + PLABEL(I,2)(6:12) (e.g.) AgCms
        if is_different_atom(ao, function_label):
            # Different atom
            ao.reset()
            ao.idx_within_same_atom = get_idx_within_the_same_atom()
            # ao was reset, so set the current subshell and gto_type again
            ao.current_ao.update(atom, subshell, gto_type)
            ao.prev_ao = copy.deepcopy(ao.current_ao)

        debug_print(f"function_label: {function_label}, ao: {ao}, idx_within_same_atom: {ao.idx_within_same_atom}")
        ao.function_types.add(function_label)

        return Function(component_func, symmetry, atom, gto_type, ao.idx_within_same_atom, num_functions, multiplicity)

    def update_symmetry_orbitals_summary(line_str: str) -> None:
        tmp_list = [0]  # Add 0 to the beginning of the list
        tmp_list.extend(int(idx) for idx in line_str[line_str.find(":") + 1 :].split())
        indices = tuple(tmp_list)
        if "large" in line_str:
            orb_summary.large_sym = indices
        elif "small" in line_str:
            orb_summary.small_sym = indices
        else:
            orb_summary.all_sym = indices

    def add_function(func: Function) -> None:
        # Create an empty dictionary if the key does not exist
        functions_info.setdefault(component_func, OrderedDict()).setdefault(symmetry, OrderedDict()).setdefault(func.atom, OrderedDict())
        # Add function information
        if func.idx_within_same_atom not in functions_info[component_func][symmetry][func.atom].keys():
            label = symmetry + func.atom
            prev_idx_per_sym = fn_summary.dirac21.get_current_func_idx(symmetry)
            prev_idx_all_sym = fn_summary.dirac19.cur_func_idx
            fn_idx_per_sym = FuncIndices(first=prev_idx_per_sym + 1, last=prev_idx_per_sym + func.num_functions)
            fn_idx_all_sym = FuncIndices(first=prev_idx_all_sym + 1, last=prev_idx_all_sym + func.num_functions)
            functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom] = AtomInfo(
                func.idx_within_same_atom, label, func.multiplicity, fn_idx_per_sym, fn_idx_all_sym
            )
        functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].add_function(func.gto_type, func.num_functions)

    def update_last_indices(func: Function) -> None:
        fn_summary.dirac19.add_num_functions(func.num_functions)
        fn_summary.dirac21.add_num_functions(symmetry, func.num_functions)
        functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].func_idx_all_sym.last = fn_summary.dirac19.cur_func_idx
        functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].func_idx_per_sym.last = fn_summary.dirac21.get_current_func_idx(symmetry)

    def check_symmetry_orbitals_summary() -> None:
        if len(orb_summary.all_sym) == 0:
            msg = "error: The number of symmetry orbitals is not found."
            raise ValueError(msg)
        elif len(orb_summary.all_sym) != len(orb_summary.large_sym) or len(orb_summary.all_sym) != len(orb_summary.small_sym):
            msg = "error: The number of elements in orb_summary.all_sym, orb_summary.large_sym, and orb_summary.small_sym must be the same."
            raise ValueError(msg)

    start_symmetry_orbitals_section = False
    component_func = "large"  # "large" or "small"
    symmetry = ""
    idx_symmetry = -1
    functions_info = FunctionsInfo()
    fn_summary = FuncNumSummary()
    orb_summary = SymmetryOrbitalsSummary()
    for line_str in dirac_output:
        words: List[str] = space_separated_parsing(line_str)
        if len(line_str) == 0:
            continue
        elif not start_symmetry_orbitals_section:
            start_symmetry_orbitals_section = is_start_symmetry_orbitals_section(words)
        elif "Number of" in line_str:
            update_symmetry_orbitals_summary(line_str)
        elif "component functions" in line_str:
            check_symmetry_orbitals_summary()
            component_func = get_component_func(line_str)
            idx_symmetry = -1
        elif "Symmetry" in line_str:
            fn_summary.dirac19.check_prev_sym_idx()
            symmetry = get_symmetry(words)
            idx_symmetry = get_symmetry_idx(line_str)
            check_symmetry_idx(idx_symmetry, symmetry, orb_summary)
            if component_func == "large":
                start = sum(orb_summary.all_sym[:idx_symmetry])
                last = sum(orb_summary.all_sym[:idx_symmetry]) + orb_summary.large_sym[idx_symmetry]
                fn_summary.dirac19.set_indices(start, last)
            elif component_func == "small":
                start = sum(orb_summary.all_sym[:idx_symmetry]) + orb_summary.large_sym[idx_symmetry]
                last = sum(orb_summary.all_sym[: idx_symmetry + 1])
                fn_summary.dirac19.set_indices(start, last)
        elif "functions" in line_str:
            func = read_func_info(words, line_str, component_func, symmetry)
            add_function(func)
            update_last_indices(func)

        elif all(char in "* \r\n" for char in line_str) and len(re.findall("[*]", line_str)) > 0:
            # all characters in line_str are * or space or line break and at least one * is included
            break  # Stop reading symmetry orbitals
    if not start_symmetry_orbitals_section:
        msg = 'ERROR: The "Symmetry Orbitals" section, which is one of the essential information sections for this program, \
is not in the DIRAC output file.\n\
Please check your DIRAC output file.\n\
Perhaps you explicitly set the .PRINT option to a negative number in one of the sections?'
        raise Exception(msg)

    return functions_info
