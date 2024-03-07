import copy
import re
from collections import OrderedDict
from io import TextIOWrapper
from typing import List, Tuple
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.atoms import AtomInfo, FuncIndices, ao, is_different_atom
from sum_dirac_dfcoef.utils import debug_print, space_separated_parsing


class LastFuncNum:
    gerade: int = 0  # Number of functions with gerade symmetry
    ungerade: int = 0  # Number of functions with ungerade symmetry
    no_inv_sym: int = 0  # Number of functions with no inversion symmetry
    sum_all: int = 0  # Sum of all functions

    def add(self, symmetry: str, num_functions: int) -> None:
        if symmetry[-1].lower() == "g":
            self.gerade += num_functions
        elif symmetry[-1].lower() == "u":
            self.ungerade += num_functions
        else:
            self.no_inv_sym += num_functions
        self.sum_all += num_functions

    def get(self, symmetry: str) -> int:
        if symmetry[-1].lower() == "g":
            return self.gerade
        elif symmetry[-1].lower() == "u":
            return self.ungerade
        else:
            return self.no_inv_sym

    def get_all(self) -> int:
        return self.sum_all


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
    #                     last_func_num: 6
    #                 }
    #            },
    #            "3": {
    #                 AtomInfo: {
    #                     mul: 2,
    #                     functions: {
    #                         "s": 3,
    #                         "p": 3,
    #                     },
    #                     last_func_num: 12
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


def get_functions_info(dirac_output: TextIOWrapper) -> FunctionsInfo:
    def is_start_symmetry_orbitals_section(words: List[str]) -> bool:
        # ref: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/dirac/dirtra.F#L3654
        if len(words) == 2 and words[0] == "Symmetry" and words[1] == "Orbitals":
            return True
        return False

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

    def check_symmetry_idx(idx_symmetry: int, symmetry: str) -> None:
        if idx_symmetry < 0:
            msg = f"error: The symmetry index must be >= 0, but got {idx_symmetry} in symmetry: {symmetry}"
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
            summary.large_sym = indices
        elif "small" in line_str:
            summary.small_sym = indices
        else:
            summary.all_sym = indices

    def check_symmetry_orbitals_summary() -> None:
        if len(summary.all_sym) == 0:
            msg = "error: The number of symmetry orbitals is not found."
            raise ValueError(msg)
        elif len(summary.all_sym) != len(summary.large_sym) or len(summary.all_sym) != len(summary.small_sym):
            msg = "error: The number of elements in summary.all_sym, summary.large_sym, and summary.small_sym must be the same."
            raise ValueError(msg)

    start_symmetry_orbitals_section = False
    component_func = ""  # "large" or "small"
    symmetry = ""
    idx_symmetry = -1
    cur_orb_idx = 0
    last_orb_idx = 0
    last_func_num = LastFuncNum()
    functions_info = FunctionsInfo()
    summary = SymmetryOrbitalsSummary()
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
            component_func = "large" if "Large" in line_str else ("small" if "Small" in line_str else "")
            idx_symmetry = -1
        elif "Symmetry" in line_str:
            if cur_orb_idx != last_orb_idx:
                msg = f"cur_orb_idx: {cur_orb_idx}, last_orb_idx: {last_orb_idx} must be the same when reading the beginning of the symmetry orbitals."
                raise ValueError(msg)
            symmetry = get_symmetry(words)
            idx_symmetry = get_symmetry_idx(line_str)
            check_symmetry_idx(idx_symmetry, symmetry)
            if component_func == "large":
                cur_orb_idx = sum(summary.all_sym[:idx_symmetry])
                last_orb_idx = sum(summary.all_sym[:idx_symmetry]) + summary.large_sym[idx_symmetry]
            elif component_func == "small":
                cur_orb_idx = sum(summary.all_sym[:idx_symmetry]) + summary.large_sym[idx_symmetry]
                last_orb_idx = sum(summary.all_sym[: idx_symmetry + 1])
        elif "functions" in line_str:
            func = read_func_info(words, line_str, component_func, symmetry)
            # Create an empty dictionary if the key does not exist
            functions_info.setdefault(component_func, OrderedDict()).setdefault(symmetry, OrderedDict()).setdefault(func.atom, OrderedDict())
            # Add function information
            if func.idx_within_same_atom not in functions_info[component_func][symmetry][func.atom].keys():
                label = symmetry + func.atom
                prev_num_per_sym = last_func_num.get(symmetry)
                ln = last_func_num.get(symmetry)
                fn_idx_per_sym = FuncIndices(first=prev_num_per_sym + 1, last=prev_num_per_sym + func.num_functions)
                fn_idx_all_sym = FuncIndices(first=cur_orb_idx + 1, last=cur_orb_idx + func.num_functions)
                functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom] = AtomInfo(
                    func.idx_within_same_atom, label, func.multiplicity, ln, fn_idx_per_sym, fn_idx_all_sym
                )
            last_func_num.add(symmetry, func.num_functions)
            cur_orb_idx += func.num_functions
            functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].add_function(func.gto_type, func.num_functions)

            functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].func_idx_all_sym.last = cur_orb_idx
            functions_info[component_func][symmetry][func.atom][func.idx_within_same_atom].func_idx_per_sym.last = last_func_num.get(symmetry)
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
