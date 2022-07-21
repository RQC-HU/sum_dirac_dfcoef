import argparse
import re
import numpy as np
import sys


class Coefficients:
    norm_const_sum: float = 0.0
    sum_of_mo_coefficient: float = 0.0
    mo_coefficient_list: "list[float]" = []

    def __repr__(self) -> str:
        return f"norm_const_sum: {self.norm_const_sum}, sum_of_mo_coefficient: {self.sum_of_mo_coefficient}, mo_coefficient_list: {self.mo_coefficient_list}"

    def reset(self):
        self.norm_const_sum = 0.0
        self.sum_of_mo_coefficient = 0.0
        self.mo_coefficient_list: "list[float]" = []


class Atoms:
    atom_nums: "list[int]" = list()
    atom_types: "list[str]" = list()
    atom_list: "list[str]" = list()
    orbital_types: "list[str]" = list()
    elements: "list[str]" = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"]

    def __repr__(self) -> str:
        return f"atom_nums: {self.atom_nums}, atom_types: {self.atom_types}, atom_list: {self.atom_list}, orbital_types: {self.orbital_types}"

    def reset(self):
        self.atom_list: "list[str]" = list()
        self.orbital_types: "list[str]" = list()


def parse_args() -> "argparse.Namespace":
    parser = argparse.ArgumentParser(description="Summarize vector priors.")
    parser.add_argument("-f", "--file", type=str, help="(required) file name of DIRAC output")
    parser.add_argument(
        "-mol",
        "--molecule",
        type=str,
        help="(required) molecule specification. Write the molecular formula (e.g. Cu2O)",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution",
    )
    return parser.parse_args()


def get_dirac_filename(args: "argparse.Namespace") -> str:
    if not args.file:
        sys.exit("ERROR: DIRAC output file is not given. Please use -f option.")
    return args.file


def get_atom_num(num_str: "str | None", elem: str) -> int:
    if num_str is None:
        return 1
    elif num_str not in elem:
        sys.exit(f"ERROR: {num_str} is not in {elem}")
    return int(num_str)


def parse_atom_num(line: str) -> int:
    num_str = "".join(filter(str.isdigit, line)) or None
    return get_atom_num(num_str, line)


def parse_atom(input_list: "list[str]", line: str, atoms: Atoms) -> str:
    letter = "".join(filter(str.isalpha, line)) or None
    if letter not in atoms.elements:
        sys.exit(f"ERROR: {line} is invalid.")
    if letter in atoms.atom_types:
        sys.exit(f"ERROR: {letter} is duplicated.\nYour input: {''.join(input_list)}")
    return letter


def parse_atom_and_the_number_of_atoms(split_list: "list[str]", atoms: Atoms) -> Atoms:
    for elem in split_list:

        letter = parse_atom(split_list, elem, atoms)

        atoms.atom_types.append(letter)

        num = parse_atom_num(elem)
        atoms.atom_nums.append(num)

        print(f"{letter} {num}")
    return atoms


def parse_molecule_input(args: "argparse.Namespace") -> Atoms:
    atoms = Atoms()
    if not args.molecule:
        sys.exit("ERROR: Molecule is not specified. Please use -mol option. (e.g. -mol Cu2O)")
    split_per_atom = re.findall("[A-Z][^A-Z]*", args.molecule)
    return parse_atom_and_the_number_of_atoms(split_per_atom, atoms)


def parse_not_empty(line: str) -> "list[str]":
    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def create_data_per_atom(
    atoms: Atoms,
    coefficients: Coefficients,
) -> "np.ndarray":
    data_type = [("atom", "U2"), ("orbital_type", "U20"), ("mo_percentage", "f8")]
    data_per_atom = np.zeros(len(atoms.orbital_types), dtype=data_type)
    for idx, (orb, var, atom) in enumerate(zip(atoms.orbital_types, coefficients.mo_coefficient_list, atoms.atom_list)):
        data_per_atom[idx]["atom"] = atom
        data_per_atom[idx]["orbital_type"] = orb
        atom_num = atoms.atom_nums[atoms.atom_types.index(atom)]
        data_per_atom[idx]["mo_percentage"] = var * 100 / (coefficients.norm_const_sum * atom_num)
    return data_per_atom


def calculate_sum_of_mo_coefficient(coefficients: Coefficients) -> float:
    return (sum([c for c in coefficients.mo_coefficient_list])) / coefficients.norm_const_sum


def write_results(
    data_per_atom: "np.ndarray",
    atoms: Atoms,
    coefficients: Coefficients,
    threshold: float,
):
    for d in data_per_atom:
        if d["mo_percentage"] > threshold:
            for _ in range(atoms.atom_nums[atoms.atom_types.index(d["atom"])]):
                orb_type = str(d["orbital_type"]).ljust(11, " ")
                print(f"{orb_type}: {d['mo_percentage']}\t%")
    print(f"Normalization constant is {coefficients.norm_const_sum}")
    print(f"sum of coefficient {coefficients.sum_of_mo_coefficient}\n")


def need_to_skip_this_line(words: "list[str]", is_read_value: bool) -> bool:
    if is_read_value and len(words) == 0:
        return False
    elif len(words) <= 1:
        return True
    else:
        return False


def need_to_write_results(words: "list[str]", is_read_value: bool) -> bool:
    if is_read_value and len(words) == 0:
        return True
    else:
        return False


def main() -> None:

    args: "argparse.Namespace" = parse_args()
    dirac_file: str = get_dirac_filename(args)
    atoms: Atoms = parse_molecule_input(args)
    threshold: float = args.threshold if args.threshold else 0.1
    is_read_value: bool = False
    symmetry_type: str = ""
    eigenvalue_num: int = 0
    eg_type: str = "E1g"
    coefficients = Coefficients()

    print(f"threshold is {threshold} %\n")
    with open(dirac_file) as f:
        for idx, line in enumerate(f):

            words: "list[str]" = parse_not_empty(line)

            if need_to_skip_this_line(words, is_read_value):
                continue

            # If this condition is true, end print eigenvalues under this line.
            # So we need to set down is_read_value flag to quit to read values.
            if need_to_write_results(words, is_read_value):
                is_read_value = False

                data_per_atom = create_data_per_atom(atoms, coefficients)

                coefficients.sum_of_mo_coefficient = calculate_sum_of_mo_coefficient(coefficients)

                write_results(
                    data_per_atom,
                    atoms,
                    coefficients,
                    threshold,
                )

                # Reset variables
                coefficients.reset()
                atoms.reset()

                continue

            # Read line and get values that we need.
            if is_read_value:

                word_len = len(words)
                atom_type = ""
                if word_len == 9:  # Ag pattern
                    symmetry_type = words[2]
                    idx = 3
                    atom_type = words[idx]  # name of atom (e.g. U)
                elif word_len == 8:
                    idx = 2
                    """
                    (e.g) words[idx] = "B3gO"
                          symmetry_type = "B3g"
                          atom_type = "O"
                    """
                    symmetry_type = words[idx][:3]
                    atom_type = words[idx][3:]

                orbital_type: str = words[idx + 1]  # orbital type (e.g. dxy)
                alpha1: float = float(words[idx + 2])
                alpha2: float = float(words[idx + 3])
                beta1: float = float(words[idx + 4])
                beta2: float = float(words[idx + 5])
                atom_orb_type = symmetry_type + atom_type + orbital_type
                if atom_orb_type not in atoms.orbital_types:
                    atoms.orbital_types.append(atom_orb_type)
                    atoms.atom_list.append(atom_type)
                    coefficients.mo_coefficient_list.append(0.0)
                coefficient = 0.0
                if atom_type not in atoms.atom_types:
                    sys.exit(f"ERROR: atom type {atom_type} is not defined. Please check your --molecule or -mol option.")
                else:
                    magnification = atoms.atom_nums[atoms.atom_types.index(atom_type)]
                    coefficient = magnification * alpha1**2 + alpha2**2 + beta1**2 + beta2**2
                coefficients.norm_const_sum += coefficient
                orb_type_idx = atoms.orbital_types.index(atom_orb_type)
                coefficients.mo_coefficient_list[orb_type_idx] += coefficient

            # If this condition is true, start print eigenvalues under this line.
            # So set up is_read_value flag to read values under this line.
            if not is_read_value and words[1] == "Electronic" and len(words) == 6:
                is_read_value = True
                if eigenvalue_num >= int(words[4][:-1]):
                    eg_type = "E1u"
                eigenvalue_num = int(words[4][:-1])
                print("Electronic no.", eigenvalue_num, eg_type)


if __name__ == "__main__":
    main()
