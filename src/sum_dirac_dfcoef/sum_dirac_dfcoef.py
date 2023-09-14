#!/usr/bin/env python3

import argparse
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


class Coefficients:
    norm_const_sum: float = 0.0
    sum_of_mo_coefficient: float = 0.0
    mo_coefficient_list: "list[float]" = list()
    orbital_types: "list[str]" = list()
    atom_list: "list[str]" = list()

    def __repr__(self) -> str:
        return f"norm_const_sum: {self.norm_const_sum}, sum_of_mo_coefficient: {self.sum_of_mo_coefficient}, mo_coefficient_list: {self.mo_coefficient_list}"

    def reset(self):
        self.norm_const_sum = 0.0
        self.sum_of_mo_coefficient = 0.0
        self.mo_coefficient_list: "list[float]" = list()
        self.orbital_types: "list[str]" = list()
        self.atom_list: "list[str]" = list()


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


class Data_per_MO:
    mo_info: str = ""
    mo_energy: float = 0.0
    data_per_orbital_types: "list[Data_per_orbital_types]" = list()
    norm_constant: float = 0.0
    sum_coefficients: float = 0.0

    def __init__(
        self,
        mo_info: str,
        mo_energy: float,
        data_per_orbital_types: "list[Data_per_orbital_types]",
        norm_constant: float,
        sum_coefficients: float,
    ) -> None:
        self.mo_info = mo_info
        self.mo_energy = mo_energy
        self.data_per_orbital_types = data_per_orbital_types
        self.norm_constant = norm_constant
        self.sum_coefficients = sum_coefficients

    def __repr__(self) -> str:
        return f"mo_info: {self.mo_info}, mo_energy: {self.mo_energy}, coefficients: {self.data_per_orbital_types}"


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
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution", dest="threshold")
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
        elif "functions" in line_str:
            # ref: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirtra.F#L3697-3699
            try:
                num_functions = int(words[0])  # ILAB(1,I)
            except (ValueError, TypeError):
                num_functions = -1  # Impossible number of functions to detect that we cannot get the number of functions from this line_str
            after_functions = line_str[line_str.find("functions:") + len("functions:") :].strip()  # PLABEL(I,2)(6:12),1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG)
            function_label = after_functions[:7].strip()  # PLABEL(I,2)(6:12)
            ao.current_subshell = function_label[3]  # e.g. "g" in "Cm g400"
            if function_label in ao.function_types or is_reverse_subshell():
                # Different atom
                ao.function_types.clear()
            multiplicity_label = after_functions[7:].strip()  # 1,(CHRSGN(NINT(CTRAN(II,K))),K,K=2,NDEG) (e.g.) 1+2+3+4
            multiplicity = len(re.findall("[+-]", multiplicity_label)) + 1  # (e.g.) 1+2=>2, 1+2+3=>3, 1+2-3-4=>4
            ao.function_types.add(function_label)
            if current_symmetry not in functions_info[current_component_function]:
                functions_info[current_component_function][current_symmetry] = dict()
            functions_info[current_component_function][current_symmetry][function_label] = {"functions": num_functions, "mul": multiplicity}
        # all characters in line_str are * or space
        elif len(re.findall("[* ]", line_str)) == len(line_str):
            break  # Stop reading symmetry orbitals

    if not start_symmetry_orbitals_section:
        raise Exception(
            "ERROR: The \"Symmetry Orbitals\" section, which is one of the essential information sections for this program,\
is not in the DIRAC output file.\n\
Please check your DIRAC output file.\n\
Perhaps you explicitly set the .PRINT option to a negative number in one of the sections?"
        )

    return functions_info


def get_coefficient(words: "list[str]", atoms: Atoms, coefficients: Coefficients, elements: "list[str]") -> None:
    """
    Nested functions to get coefficient

    words is a list of strings that is split by space.
    words[0] is the number of MO.
    words[1] is the L or S, Large or Small.
    words[2:-5] are symmetry type, Atom and Orbital type.
      Sometimes these elements cannot be separated because there is no space available <--- This is the reason why we use words[2:-5].
    words[-4:] is the coefficient.

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
    def get_types() -> "tuple[str, str, str]":
        sym_and_atom_and_orb_str = " ".join(words[2:-4])
        splitted_by_capital = re.findall("[A-Z][^A-Z]*", sym_and_atom_and_orb_str)
        symmetry_type = splitted_by_capital[0].strip()
        # atom_type must be 1 or 2 letters.
        # check splitted_by_capital[1][:2] is not out of range
        if len(splitted_by_capital[1]) >= 2 and splitted_by_capital[1][:2] in elements:  # 2 letters (e.g. Cu)
            atom_type = splitted_by_capital[1][:2]
            orbital_type = splitted_by_capital[1][2:]
        elif splitted_by_capital[1][0] in elements:  # 1 letter (e.g. C)
            atom_type = splitted_by_capital[1][0]
            orbital_type = splitted_by_capital[1][1:]
        else:
            sys.exit(f"ERROR: {splitted_by_capital[1][:1]} is invalid atom type.")

        # orbital_type does not have space or numbers.
        orbital_type = orbital_type.lstrip("0123456789 ")

        # Return symmetry_type, atom_type, orbital_type with no space.
        return symmetry_type.strip(), atom_type.strip(), orbital_type.strip()

    def add_orbital_type(atom_type: str, atom_orb_type: str) -> None:
        if atom_orb_type not in coefficients.orbital_types:
            coefficients.orbital_types.append(atom_orb_type)
            coefficients.atom_list.append(atom_type)
            coefficients.mo_coefficient_list.append(0.0)

    def isfloat(parameter):
        if not parameter.isdecimal():
            try:
                float(parameter)
                return True
            except ValueError:
                return False
        else:
            return False

    def get_coefficient() -> float:
        """
        (e.g)
        words = ["g400", "0.0000278056", "0.0000000000", "0.0000000000", "0.0000000000"]
        """
        alpha1: float = float(words[-4]) if isfloat(words[-4]) else 0.0
        alpha2: float = float(words[-3]) if isfloat(words[-3]) else 0.0
        beta1: float = float(words[-2]) if isfloat(words[-2]) else 0.0
        beta2: float = float(words[-1]) if isfloat(words[-1]) else 0.0
        return alpha1**2 + alpha2**2 + beta1**2 + beta2**2

    def add_coefficient(coefficient: float, atom_orb_type: str) -> None:
        magnification = atoms.atom_nums[atoms.atom_types.index(atom_type)]

        coefficient = magnification * coefficient
        coefficients.norm_const_sum += coefficient

        orb_type_idx = coefficients.orbital_types.index(atom_orb_type)
        coefficients.mo_coefficient_list[orb_type_idx] += coefficient

    def parse_words(words: "list[str]") -> "list[str]":
        new_words = []
        # Parses multiple numbers that are sometimes connected without separating them with spaces
        # In DIRAC version <22.0 the number of decimal places is fixed at 10 (decimal_num = 10)
        for word in words:
            num_of_dots = word.count(".")
            if num_of_dots == 0 or num_of_dots == 1:
                new_words.append(word)
            elif num_of_dots >= 2:
                decimal_num = 10
                dotidx = word.find(".")
                while dotidx != -1:
                    word2 = word[: dotidx + decimal_num + 1]
                    new_words.append(word2)
                    word = word[dotidx + decimal_num + 1 :]
                    dotidx = word.find(".")

        return new_words

    """
    Main function to get coefficient
    """
    words = parse_words(words)
    symmetry_type, atom_type, orbital_type = get_types()
    coefficient = get_coefficient()
    atom_orb_type = symmetry_type + atom_type + orbital_type

    add_orbital_type(atom_type, atom_orb_type)
    add_coefficient(coefficient, atom_orb_type)

    return None


def create_results_for_current_mo(args: "argparse.Namespace", atoms: Atoms, coefficients: Coefficients) -> "tuple[list[Data_per_orbital_types], float, float]":
    """
    Nested functions to create results for current MO
    """

    def create_data_per_orbital_types():
        data_per_orbital_types: "list[Data_per_orbital_types]" = []
        for orb, coefficient, atom in zip(
            coefficients.orbital_types,
            coefficients.mo_coefficient_list,
            coefficients.atom_list,
        ):
            atom_num = atoms.atom_nums[atoms.atom_types.index(atom)]
            data = Data_per_orbital_types(
                atom=atom,
                orbital_type=orb,
                mo_percentage=coefficient * 100 / (coefficients.norm_const_sum * atom_num),
            )
            if data.mo_percentage >= args.threshold:
                for _ in range(atom_num):
                    data_per_orbital_types.append(data)
        return data_per_orbital_types

    def calculate_sum_of_mo_coefficient() -> float:
        return (sum([c for c in coefficients.mo_coefficient_list])) / coefficients.norm_const_sum

    """
    Main function to create results for current MO
    """
    data_per_orbital_types = create_data_per_orbital_types()
    data_per_orbital_types.sort(key=lambda x: x.mo_percentage, reverse=True)
    normalization_constant = 0.0
    sum_of_coefficient = 0.0
    if args.debug:
        normalization_constant = coefficients.norm_const_sum
        sum_of_coefficient = calculate_sum_of_mo_coefficient()

    return data_per_orbital_types, normalization_constant, sum_of_coefficient


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
        output_name = f"{args.mol}.out"
        output_path = os.path.join(os.getcwd(), output_name)
    else:
        output_name = args.output
        output_path = os.path.abspath(output_name)
    return output_path


def write_results(args: "argparse.Namespace", file: TextIOWrapper, data_all_mo: "list[Data_per_MO]") -> None:
    """
    Write results to stdout
    """

    for mo in data_all_mo:
        digit_int = len(str(int(mo.mo_energy)))  # number of digits of integer part
        # File write but if args.compress is True \n is not added
        mo_info_energy = f"{mo.mo_info} {mo.mo_energy:{digit_int}.{args.decimal}f}" + ("\n" if not args.compress else "")
        file.write(mo_info_energy)

        d: Data_per_orbital_types
        for d in mo.data_per_orbital_types:
            if args.compress:
                orb_type = str(d.orbital_type)
                output_str = f" {orb_type} {d.mo_percentage:.{args.decimal}f}"
                file.write(output_str)
            else:
                orb_type = str(d.orbital_type).ljust(11, " ")
                output_str = f"{orb_type} {d.mo_percentage:{args.decimal+4}.{args.decimal}f} %\n"
                file.write(output_str)
        file.write("\n")  # add empty line
        debug_print_wrapper(args, f"Normalization constant is {mo.norm_constant:.{args.decimal}f}")
        debug_print_wrapper(args, f"sum of coefficient {mo.sum_coefficients:.{args.decimal}f}")


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
    # fmt: off
    elements: "list[str]" = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    # fmt: on

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
                    (
                        data,
                        norm_constant,
                        sum_coefficients,
                    ) = create_results_for_current_mo(args, atoms, coefficients)
                    if args.compress:
                        info = f"{mo_sym_type} {electron_number}"
                    else:
                        info = f"Electronic no. {electron_number} {mo_sym_type}"
                    if is_electronic:
                        data_all_electronic_mo.append(
                            Data_per_MO(
                                mo_info=info,
                                mo_energy=mo_energy,
                                data_per_orbital_types=data,
                                norm_constant=norm_constant,
                                sum_coefficients=sum_coefficients,
                            )
                        )
                    else:  # Positronic
                        data_all_positronic_mo.append(
                            Data_per_MO(
                                mo_info=info,
                                mo_energy=mo_energy,
                                data_per_orbital_types=data,
                                norm_constant=norm_constant,
                                sum_coefficients=sum_coefficients,
                            )
                        )
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
                except ValueError:  # If *** is printed, we have no information about what number this MO is. Therefore, we assume that electron_number is the next number after prev_electron_number.
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
                get_coefficient(words, atoms, coefficients, elements)
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
