#!/usr/bin/env python3

import argparse
import bisect
import copy
from io import TextIOWrapper
import os
import re
import sys
from typing import Annotated

from annotated_types import MaxLen
from pydantic import BaseModel


class Atoms:
    total_number_of_atoms: int = 0
    number_of_atom_types: int = 0
    atom_info: "dict[str, int]" = dict()  # key: atom type, value: number of atoms

    def __init__(self, atom_info: "dict[str, int]") -> None:
        self.total_number_of_atoms = sum(atom_info.values())
        self.number_of_atom_types = len(atom_info)
        self.atom_info = atom_info

    def __repr__(self) -> str:
        return f"total_number_of_atoms: {self.total_number_of_atoms}, number_of_atom_types: {self.number_of_atom_types}, atom_info: {self.atom_info}"


class AtomicOrbitals(BaseModel, validate_assignment=True):
    prev_subshell: Annotated[str, MaxLen(max_length=1)] = "s"
    current_subshell: Annotated[str, MaxLen(max_length=1)] = "s"
    function_types: "set[str]" = set()

    def reset(self):
        self.prev_subshell = "s"
        self.current_subshell = "s"
        self.function_types.clear()


class AOFunction(BaseModel, validate_assignment=True):
    atom_type: str
    orbital_type: str
    symmetry_type: str
    function_num: int
    mul: int


class SymmetryOrbital(BaseModel, validate_assignment=True):
    # {"Ag": dict(), "B1g": dict(), ...}
    name: str
    function_labels: "dict[str, AOFunction]" = dict()


class FunctionInfo(BaseModel, validate_assignment=True):
    # {"large orbitals": SymmetryOrbital(), "small orbitals": SymmetryOrbital()}
    name: str
    symmetry_orbitals: SymmetryOrbital


class Coefficient(BaseModel, validate_assignment=True):
    vector_num: int
    function_label: str
    # atom: str
    coefficient: float
    magnification: int

    def __repr__(self) -> str:
        super().__repr__()
        return f"vector_num: {self.vector_num}, function_label: {self.function_label}, coefficient: {self.coefficient}, magnification: {self.magnification}"


class Coefficients:
    norm_const_sum: float = 0.0
    coef_dict: "dict[str, Coefficient]" = dict()
    coef_list: "list[Coefficient]" = list()

    def __repr__(self) -> str:
        # return f"norm_const_sum: {self.norm_const_sum}, coef_list: {self.coef_list}"
        return f"norm_const_sum: {self.norm_const_sum}, coef: {[coef.coefficient for coef in self.coef_list]}"

    def add_coefficient(self, coef: Coefficient) -> None:
        # self.coef_list.append(coef)
        if coef.function_label in self.coef_dict:
            self.coef_dict[coef.function_label].coefficient += coef.coefficient
        else:
            self.coef_dict[coef.function_label] = coef
        self.norm_const_sum += coef.coefficient * coef.magnification

    def reset(self):
        self.norm_const_sum = 0.0
        self.coef_dict.clear()
        self.coef_list.clear()


class Data_per_MO:
    coefficients: Coefficients
    mo_energy: float
    mo_info: str

    def __init__(self, coefficients: Coefficients, mo_energy: float, mo_info: str) -> None:
        self.coefficients = coefficients
        self.mo_energy = mo_energy
        self.mo_info = mo_info


class VectorInfo:
    vector_num: int
    atom_label: str
    function_label: str
    subshell: str
    coefficient: float


class Data_per_orbital_types:
    atom: str = ""
    orbital_type: str = ""
    mo_percentage: float = 0.0

    def __init__(self, atom: str, orbital_type: str, mo_percentage: float) -> None:
        self.atom = atom
        self.orbital_type = orbital_type
        self.mo_percentage = mo_percentage

    def __repr__(self) -> str:
        return f"atom: {self.atom}, orbital_type: {self.orbital_type}, mo_percentage: {self.mo_percentage}"

    def reset(self):
        self.atom = ""
        self.orbital_type = ""
        self.mo_percentage = 0.0


class PrintVersionExitAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        from .__about__ import __version__

        print(f"{__version__}")
        exit()


def parse_args() -> "argparse.Namespace":
    parser = argparse.ArgumentParser(
        description="Summarize the coefficients from DIRAC output file that *PRIVEC option is used. (c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)"
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="(required) file name of DIRAC output", dest="file")
    # parser.add_argument("-m", "--mol", type=str, required=True, help="(required) molecule specification. Write the molecular formula (e.g. Cu2O). ** DON'T write the rational formula (e.g. CH3OH) **")
    parser.add_argument("-o", "--output", type=str, help="Output file name. Default: (-m or --mol option value).out (e.g) --m H2O => print to H2O.out", dest="output")
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        help="Compress output. Display all coefficients on one line for each MO. This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.",
        dest="compress",
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.1, help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution", dest="threshold"
    )
    parser.add_argument(
        "-d",
        "--decimal",
        type=int,
        default=5,
        choices=range(1, 16),
        help="Set the decimal places. Default: 5 (e.g) --decimal=3 => print orbital with 3 decimal places (0.123, 2.456, ...). range: 1-15",
        dest="decimal",
    )
    parser.add_argument("-a", "--all-write", action="store_true", help="Print all MOs(Positronic and Electronic).", dest="all_write")
    parser.add_argument("-p", "--positronic-write", action="store_true", help="Print only Positronic MOs.", dest="positronic_write")
    parser.add_argument("-v", "--version", action=PrintVersionExitAction, help="Print version and exit", dest="version")
    parser.add_argument("--debug", action="store_true", help="print debug output (Normalization constant, Sum of MO coefficient)", dest="debug")
    parser.add_argument("--no-sort", action="store_true", help="Don't sort the output by MO energy")
    # If -v or --version option is used, print version and exit
    return parser.parse_args()


def debug_print_wrapper(args: "argparse.Namespace", str: str):
    # print debug message if --debug option is used
    if args.debug:
        print(str)


def is_this_row_for_coefficients(words: "list[str]") -> bool:
    # min: 4 coefficients and other words => 5 words
    if 5 <= len(words) <= 9 and words[0].isdigit():
        return True
    else:
        return False


def need_to_skip_this_line(words: "list[str]") -> bool:
    if len(words) <= 1:
        return True
    else:
        return False


def need_to_create_results_for_current_mo(words: "list[str]", is_reading_coefficients: bool) -> bool:
    if is_reading_coefficients and len(words) <= 1:
        return True
    else:
        return False


def need_to_get_mo_sym_type(words: "list[str]", start_mo_coefficients: bool) -> bool:
    if not start_mo_coefficients and len(words) == 3 and words[0] == "Fermion" and words[1] == "ircop":
        return True
    return False


def need_to_start_mo_section(words: "list[str]", start_mo_coefficients: bool) -> bool:
    if not start_mo_coefficients and words[1] == "Electronic" and words[2] == "eigenvalue" and "no." in words[3]:
        return True
    elif not start_mo_coefficients and words[1] == "Positronic" and words[2] == "eigenvalue" and "no." in words[3]:
        return True
    return False


def get_dirac_filename(args: "argparse.Namespace") -> str:
    if not args.file:
        sys.exit("ERROR: DIRAC output file is not given. Please use -f option.")
    return args.file


def space_separated_parsing(line: str) -> "list[str]":
    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def get_atoms_and_basis_sets(dirac_output: TextIOWrapper) -> Atoms:
    """
    (e.g.)
    Atoms and basis sets
    --------------------
    Number of atom types :    5
    Total number of atoms:   57
    label    atoms   charge   prim    cont     basis
    ----------------------------------------------------------------------
    Cm          1      96     415     415      L  - [33s29p20d13f3g|33s29p20d13f3g]
    O4          4       8      34      34      L  - [10s6p1d|10s6p1d]
    N4          4       7      34      34      L  - [10s6p1d|10s6p1d]
    C          24       6      34      34      L  - [10s6p1d|10s6p1d]
    H          24       1       9       9      L  - [6s1p|6s1p]
    ----------------------------------------------------------------------
                            1719    1719      L  - large components
    ----------------------------------------------------------------------
    total:     57     324    1719    1719
    """

    def is_start_atoms_and_basis_sets(words: "list[str]") -> bool:
        # ref: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/abacus/herrdn.F#L2402
        if len(words) >= 4 and words[0] == "Atoms" and words[1] == "and" and words[2] == "basis" and words[3] == "sets":
            return True
        return False

    def is_end_atoms_and_basis_sets(words: "list[str]") -> bool:
        # ref: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/abacus/herrdn.F#L2503-2504
        if len(words) == 5 and words[0] == "total:" and words[1].isdigit() and words[2].isdigit() and words[3].isdigit() and words[4].isdigit():
            return True
        return False

    def is_include_label_and_atoms(line: "list[str]") -> bool:
        try:
            # Test for label and atoms
            # Expected format: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/abacus/herrdn.F#L2455-2457
            _ = str(line[0])  # label
            _ = int(line[1])  # atoms
            _ = int(line[2])  # charge
            _ = int(line[3])  # prim
            _ = int(line[4])  # cont
            return True
        except (ValueError, TypeError, IndexError):
            # This line does not contain label and atoms
            # Expected exceptions
            return False

    def validate_atoms_and_basis_sets_data() -> None:
        if number_of_atom_types != 0 and len(atom_info) != number_of_atom_types:
            sys.exit(
                "ERROR: Number of atom types is not equal to the number of atoms in the molecule specification.\n\
    Please check your Atoms and basis sets section in DIRAC output file.\n\
    Expected format: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/abacus/herrdn.F#L2455-2457"
            )
        if total_number_of_atoms != 0 and sum(atom_info.values()) != total_number_of_atoms:
            sys.exit(
                "ERROR: Total number of atoms is not equal to the number of atoms in the molecule specification.\n\
    Please check your Atoms and basis sets section in DIRAC output file.\n\
    Expected format: https://gitlab.com/dirac/dirac/-/blob/de590d17dd38da238ff417b4938d69564158cd7f/src/abacus/herrdn.F#L2455-2457"
            )

    number_of_atom_types = 0
    total_number_of_atoms = 0
    atom_info: "dict[str, int]" = dict()
    start_atoms_and_basis_sets = False
    for line_str in dirac_output:
        words: "list[str]" = space_separated_parsing(line_str)

        if not start_atoms_and_basis_sets:
            start_atoms_and_basis_sets = is_start_atoms_and_basis_sets(words)
            continue  # Skip this line_str and start reading labels and atoms
        # Check if this line_str is the end of Atoms and basis sets section
        elif is_end_atoms_and_basis_sets(words):
            break  # Stop reading labels and atoms

        if is_include_label_and_atoms(words):
            label = words[0]  # label
            atom_num = int(words[1])  # atoms
            atom_info[label] = atom_num
        if "Number of atom types" in line_str:
            number_of_atom_types = int(words[-1])
        elif "Total number of atoms" in line_str:
            total_number_of_atoms = int(words[-1])

    validate_atoms_and_basis_sets_data()
    return Atoms(atom_info)


def get_symmetry_orbitals(dirac_output: TextIOWrapper) -> "dict[str, dict[str, dict[str, dict[str, int]]]]":
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
    current_component_function = ""  # "large orbitals" or "small orbitals"
    # functions_info = {"large orbitals": {"Ag": {"labels: {"C  s": {"functions": 3, "mul": 2}, "C  p": {"functions": 3, "mul": 2}, ...}, "check": False}, "B1g": {...}, ...}, "small orbitals": {...}}
    functions_info: "dict[str, dict[str, dict[str, dict[str, int]]]]" = {"large orbitals": dict(), "small orbitals": dict()}
    current_symmetry = ""
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
            current_component_function = "large orbitals" if "Large" in line_str else ("small orbitals" if "Small" in line_str else "")
        elif "Symmetry" in line_str:
            current_symmetry = words[1]
            bra_idx = current_symmetry.find("(")
            if bra_idx != -1:
                current_symmetry = current_symmetry[:bra_idx]
        elif "functions" in line_str:
            # ref: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
            try:
                num_functions = int(words[0])  # ILAB(1,I)
            except (ValueError, TypeError):
                num_functions = -1  # Impossible number of functions to detect that we cannot get the number of functions from this line_str
            after_functions = line_str[line_str.find("functions:") + len("functions:") :].strip()  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)
            function_label = current_symmetry + after_functions[:7].strip()  # current_symmetry + PLABEL(I,2)(6:12)
            ao.current_subshell = after_functions[3]  # e.g. "g" in "Cm g400"
            # remove space from function_label
            function_label = function_label.replace(" ", "")
            if function_label in ao.function_types or is_reverse_subshell():
                # Different atom
                ao.function_types.clear()
            multiplicity_label = after_functions[7:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
            multiplicity = len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4
            ao.function_types.add(function_label)
            if current_symmetry not in functions_info[current_component_function]:
                functions_info[current_component_function][current_symmetry] = dict()
            functions_info[current_component_function][current_symmetry][function_label] = {"functions": num_functions, "mul": multiplicity}
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

    return functions_info


def get_coefficient(line: str, orbitals: "dict[str, dict[str, dict[str, dict[str, int]]]]") -> Coefficient:
    """
    Nested functions to get coefficient
    (e.g.)
    sym_and_atom_and_orb_str = "B3gCldyz"
    symmetry_type = "B3g"
    atom_type = "Cl"
    orbital_type = "dyz"
    """
    # ref (print ): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L388-389
    # ref (format): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L453
    # FORMAT(3X,I5,2X,A12,2X,4F14.10)
    # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L168-169

    def is_float(parameter: str):
        if not parameter.isdecimal():
            try:
                float(parameter)
                return True
            except ValueError:
                return False
        else:
            return False

    def parse_line(line: "str") -> Coefficient:
        """
        This function parses the line that contains the coefficient.

        line write source: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L440-442
        line format source: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L453
        line format : FORMAT(3X,I5,2X,A12,2X,4F14.10)
        """
        words = line.split()
        # JS (I5): Serial number of the vector
        vec_num = int(words[0])

        #
        # PLABEL(IPLAB(IBAS(IFRP)+JS,2),2) (A12): Information about the vector to identify the vector
        #
        # PLABEL source: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L168-169
        # PLABEL(NLAB,2) = CLS(IC)//' '//REP(IRP)//NAMN(ICENT)(1:3)//GTOTYP(ITYP)
        # CLS (A1): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L45
        # REP (A3): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/include/pgroup.h#L16
        # NAMN (A3, defined as A4, but only (1:3) is used): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/include/nuclei.h#L25
        # GTOTYP (A4): https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/include/ccom.h#L8
        current_component_function = "large orbitals" if line[10] == "L" else ("small orbitals" if line[10] == "S" else "")  # CLS
        symmetry_label = line[12:15].strip()  # REP (e.g. "Ag ")
        function_label = line[12:22].strip().replace(" ", "")  # REP + NAMN + GTOTYP (e.g. "Ag Cm s   " => "AgCms")

        # COEF (4F14.10)
        # coefficients = [line[24:38], line[38:52], line[52:66], line[66:80]]
        coef_num = 4
        coef_len = 14
        coef_start_idx = 24
        coefficient = sum(
            [
                pow(float(line[i : i + coef_len]), 2) if is_float(line[i : i + coef_len]) else pow(-100, 2)
                for i in range(coef_start_idx, coef_start_idx + coef_len * coef_num, coef_len)
            ]
        )

        magnification = orbitals[current_component_function][symmetry_label][function_label]["mul"]

        return Coefficient(vector_num=vec_num, function_label=function_label, coefficient=coefficient, magnification=magnification)

    """
    Main function to get coefficient
    """
    coef: Coefficient = parse_line(line)

    return coef


def check_start_vector_print(words: "list[str]") -> bool:
    # ****************************** Vector print ******************************
    if len(words) < 4:
        return False
    elif words[1] == "Vector" and words[2] == "print":
        return True
    return False


def check_end_vector_print(
    words: "list[str]",
    start_vector_print: bool,
    start_mo_section: bool,
    start_mo_coefficients: bool,
    is_reading_coefficients: bool,
) -> bool:
    # https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/issues/7#issuecomment-1377969626
    if len(words) >= 2 and start_vector_print and start_mo_section and not start_mo_coefficients and not is_reading_coefficients:
        return True
    return False


def get_output_path(args: "argparse.Namespace") -> str:
    if args.output is None:
        output_name = "sum_dirac_dfcoef.out"
        output_path = os.path.join(os.getcwd(), output_name)
    else:
        output_name = args.output
        output_path = os.path.abspath(output_name)
    return output_path


def create_results_for_current_mo(
    args: "argparse.Namespace", coefficients: Coefficients, electron_number: int, mo_energy: float, mo_sym_type: str, is_electronic: bool
) -> "Data_per_MO":
    """
    Create results for current MO
    """

    def coef_sort_and_remove(args: "argparse.Namespace", coefficients: Coefficients) -> Coefficients:
        """
        Sort coefficients and remove elements below threshold
        """

        def remove_elements_below_threshold(args: "argparse.Namespace", coefficients: Coefficients) -> Coefficients:
            """
            Remove elements below threshold
            """
            idx = bisect.bisect_left([c.coefficient / coefficients.norm_const_sum for c in coefficients.coef_list], args.threshold)
            coefficients.coef_list = coefficients.coef_list[idx:] if idx < len(coefficients.coef_list) else []
            return coefficients

        coefficients.coef_list = [v for v in coefficients.coef_dict.values()]
        if not args.no_sort:
            coefficients.coef_list.sort(key=lambda x: x.coefficient)  # ascending order
        coefficients = remove_elements_below_threshold(args, coefficients)
        coefficients.coef_list.reverse()  # descending order
        return coefficients

    coefficients = coef_sort_and_remove(args, coefficients)
    info: str
    if is_electronic:
        info = f"{mo_sym_type} {electron_number}" if args.compress else f"Electronic no. {electron_number} {mo_sym_type}"
    else:  # Positronic
        info = f"{mo_sym_type} {electron_number}" if args.compress else f"Positronic no. {electron_number} {mo_sym_type}"

    return Data_per_MO(mo_info=info, mo_energy=mo_energy, coefficients=copy.deepcopy(coefficients))


def write_results(args: "argparse.Namespace", file: TextIOWrapper, data_all_mo: "list[Data_per_MO]") -> None:
    """
    Write results to stdout
    """

    for mo in data_all_mo:
        digit_int = len(str(int(mo.mo_energy)))  # number of digits of integer part
        # File write but if args.compress is True \n is not added
        mo_info_energy = f"{mo.mo_info} {mo.mo_energy:{digit_int}.{args.decimal}f}" + ("\n" if not args.compress else "")
        file.write(mo_info_energy)

        if not args.no_sort:
            mo.coefficients.coef_list.sort(key=lambda x: x.coefficient, reverse=True)
        # print(f"mo.coefficients.coef_list: {mo.coefficients.coef_list}")
        for c in mo.coefficients.coef_list:
            percentage = c.coefficient / mo.coefficients.norm_const_sum * 100
            output_str: str
            if args.compress:
                output_str = f" {c.function_label} {percentage:.{args.decimal}f}"
            else:
                output_str = f"{c.function_label} {percentage:{args.decimal+4}.{args.decimal}f} %\n"
            for idx in range(c.magnification):
                # print(f"{idx}, {output_str}")
                file.write(output_str)
        file.write("\n")  # add empty line
        debug_print_wrapper(args, f"sum of coefficient {mo.coefficients.norm_const_sum:.{args.decimal}f}")


def main() -> None:
    start_mo_coefficients: bool = False
    start_mo_section: bool = False
    is_reading_coefficients: bool = False
    start_vector_print: bool = False
    is_electronic: bool = False
    electron_number: int = 0
    prev_electron_number: int = electron_number
    mo_energy: float = 0.0
    mo_sym_type: str = ""
    coefficients: Coefficients = Coefficients()

    args: "argparse.Namespace" = parse_args()
    dirac_filename: str = get_dirac_filename(args)
    dirac_output = open(dirac_filename, encoding="utf-8")
    atoms: Atoms = get_atoms_and_basis_sets(dirac_output)
    print(atoms)
    orbitals = get_symmetry_orbitals(dirac_output)
    data_all_electronic_mo: "list[Data_per_MO]" = []
    data_all_positronic_mo: "list[Data_per_MO]" = []
    with open(dirac_filename, encoding="utf-8") as f:
        for line in f:
            words: "list[str]" = space_separated_parsing(line)

            if not start_vector_print:
                if check_start_vector_print(words):
                    start_vector_print = True
                continue

            if need_to_get_mo_sym_type(words, start_mo_coefficients):
                mo_sym_type = words[2]

            elif need_to_skip_this_line(words):
                # if atoms is unbound, raise exception
                if need_to_create_results_for_current_mo(words, is_reading_coefficients):
                    start_mo_coefficients = False
                    data = create_results_for_current_mo(args, coefficients, electron_number, mo_energy, mo_sym_type, is_electronic)
                    if is_electronic:
                        data_all_electronic_mo.append(data)
                    else:  # Positronic
                        data_all_positronic_mo.append(data)
                    # Reset variables
                    coefficients.reset()
                    debug_print_wrapper(args, f"End of reading {electron_number}th MO")
                    is_reading_coefficients = False
                continue

            elif need_to_start_mo_section(words, start_mo_coefficients):
                """
                (e.g.)
                words = ["*", "Electronic", "eigenvalue", "no.", "22:", "-2.8417809384721"]
                words = ["*", "Electronic", "eigenvalue", "no.122:", "-2.8417809384721"]
                """
                start_mo_section = True
                start_mo_coefficients = True
                if words[1] == "Positronic":
                    is_electronic = False
                elif words[1] == "Electronic":
                    is_electronic = True
                else:
                    raise Exception("Unknown MO type")
                try:
                    electron_number = int(words[-2][:-1].replace("no.", ""))
                except ValueError:
                    # If *** is printed, we have no information about what number this MO is.
                    # Therefore, we assume that electron_number is the next number after prev_electron_number.
                    electron_number = prev_electron_number + 1
                prev_electron_number = electron_number
                mo_energy = float(words[-1])
                continue

            elif check_end_vector_print(
                words,
                start_vector_print,
                start_mo_section,
                start_mo_coefficients,
                is_reading_coefficients,
            ):
                break

            # Read coefficients or the end of coefficients section
            elif start_mo_coefficients:
                if not is_this_row_for_coefficients(words):
                    continue
                is_reading_coefficients = True
                coef = get_coefficient(line, orbitals)
                coefficients.add_coefficient(coef)
    # End of reading file
    output_path = get_output_path(args)
    file = open(output_path, "w")
    if args.all_write or args.positronic_write:
        if not args.no_sort:
            data_all_positronic_mo.sort(key=lambda x: x.mo_energy)
        write_results(args, file, data_all_positronic_mo)
        file.write("\n")  # Add a blank line
    if args.all_write or not args.positronic_write:  # Electronic
        if not args.no_sort:
            data_all_electronic_mo.sort(key=lambda x: x.mo_energy)
        write_results(args, file, data_all_electronic_mo)
    file.close()
