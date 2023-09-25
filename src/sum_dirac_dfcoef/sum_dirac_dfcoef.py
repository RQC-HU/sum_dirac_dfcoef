#!/usr/bin/env python3

import argparse
import bisect
import copy
from io import TextIOWrapper
import os
import sys

from pydantic import BaseModel

from .args import parse_args
from .atoms_and_basis_sets import get_atoms_and_basis_sets
from .utils import debug_print, space_separated_parsing
from .functions_info import FunctionsInfo, get_functions_info


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


def get_coefficient(line: str, orbitals: FunctionsInfo) -> Coefficient:
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
        component_func = "large orbitals" if line[10] == "L" else ("small orbitals" if line[10] == "S" else "")  # CLS
        symmetry_label = line[12:15].strip()  # REP (e.g. "Ag ")
        atom_label = line[15:18].strip()  # NAMN (e.g. "Cm "), atom_labe="Cm"
        gto_type = line[18:22].strip()  # GTOTYP (e.g. "s   "), gto_type="s"
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

        magnification = int(orbitals[component_func][symmetry_label][atom_label][gto_type][1]["mul"])

        return Coefficient(vector_num=vec_num, function_label=function_label, coefficient=coefficient, magnification=magnification)

    """
    Main function to get coefficient
    """
    coef = parse_line(line)

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
        debug_print(args, f"sum of coefficient {mo.coefficients.norm_const_sum:.{args.decimal}f}")


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
    atoms = get_atoms_and_basis_sets(dirac_output)
    print(atoms)
    orbitals = get_functions_info(dirac_output)
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
                    debug_print(args, f"End of reading {electron_number}th MO")
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
