#!/usr/bin/env python3

import argparse
from io import TextIOWrapper
import os
import re
import sys

from pymolecule_parser import parse


class Atoms:
    atom_nums: "list[int]" = list()
    atom_types: "list[str]" = list()

    def __init__(self, atom_nums: "list[int]", atom_types: "list[str]") -> None:
        self.atom_nums = atom_nums
        self.atom_types = atom_types

    def __repr__(self) -> str:
        return f"atom_nums: {self.atom_nums}, atom_types: {self.atom_types}"


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
    parser = argparse.ArgumentParser(description="Summarize the coefficients from DIRAC output file that *PRIVEC option is used. (c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)")
    parser.add_argument("-i", "--input", type=str, required=True, help="(required) file name of DIRAC output", dest="file")
    parser.add_argument("-m", "--mol", type=str, required=True, help="(required) molecule specification. Write the molecular formula (e.g. Cu2O). ** DON'T write the rational formula (e.g. CH3OH) **")
    parser.add_argument("-o", "--output", type=str, help="Output file name. Default: (-m or --mol option value).out (e.g) --m H2O => print to H2O.out", dest="output")
    parser.add_argument("-c", "--compress", action="store_true", help="Compress output. Display all coefficients on one line for each MO. This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.", dest="compress")
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution", dest="threshold")
    parser.add_argument("-d", "--decimal", type=int, default=5, choices=range(1, 16), help="Set the decimal places. Default: 5 (e.g) --decimal=3 => print orbital with 3 decimal places (0.123, 2.456, ...). range: 1-15", dest="decimal")
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


def parse_molecule_input(args: "argparse.Namespace") -> Atoms:
    """
    Parse the molecule input and return the Atoms object.

    """

    molecule_dict = parse(args.mol)
    atom_nums = list(molecule_dict.values())  # list of counts of atoms
    atom_types = list(molecule_dict.keys())  # list of atom types
    atoms = Atoms(atom_nums, atom_types)
    return atoms


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

    def check_atom_type(atom_type: str) -> None:
        if atom_type not in atoms.atom_types:
            print(
                "WARNING: ",
                atom_type,
                " is not in the molecule specification. This orbital will be ignored.",
            )
            print(" ".join(atoms.atom_types))
            sys.exit(f"ERROR: atom type {atom_type} is not defined. Please check your -m or --mol option.")

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
    check_atom_type(atom_type)
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
    mo_energy: float = 0.0
    mo_sym_type: str = ""
    coefficients: Coefficients = Coefficients()
    # fmt: off
    elements: "list[str]" = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    # fmt: on

    args: "argparse.Namespace" = parse_args()
    dirac_file: str = get_dirac_filename(args)
    atoms: Atoms = parse_molecule_input(args)

    data_all_electronic_mo: "list[Data_per_MO]" = []
    data_all_positronic_mo: "list[Data_per_MO]" = []
    with open(dirac_file, encoding="utf-8") as f:
        for line in f:
            words: "list[str]" = space_separated_parsing(line)
            if not start_vector_print:
                if check_start_vector_print(words):
                    start_vector_print = True
                continue

            if need_to_get_mo_sym_type(words, start_mo_coefficients):
                mo_sym_type = words[2]

            elif need_to_skip_this_line(words):
                # End of reading the specific MO coefficients
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
                electron_number = int(words[-2][:-1].replace("no.", ""))
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
