import argparse
import re
import sys


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


class Data_per_orbital_types:
    atom: "list[str]" = list()  # (e.g.) ["U","O","O",....,"U"]
    orbital_type: "list[str]" = list()  # (e.g.) ["B1uUfxxz","B1uOs","B1uOpz",....,"B2uUfyzz"]
    mo_percentage: "list[float]" = list()  # (e.g.) [0.1576884676852207,81.58228712019573,"O",....,"U"]

    def __repr__(self) -> str:
        return f"atom: {self.atom}, orbital_type: {self.orbital_type}, mo_percentage: {self.mo_percentage}"

    def reset(self):
        self.atom: "list[str]" = list()
        self.orbital_type: "list[str]" = list()
        self.mo_percentage: "list[float]" = list()


def parse_args() -> "argparse.Namespace":
    parser = argparse.ArgumentParser(description="Summarize the coefficients from DIRAC output file that *PRIVEC option is used. (c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)")
    parser.add_argument("-f", "--file", type=str, help="(required) file name of DIRAC output")
    parser.add_argument("-mol", "--molecule", type=str, help="(required) molecule specification. Write the molecular formula (e.g. Cu2O)")
    parser.add_argument("-t", "--threshold", type=float, help="threshold. Default: 0.1 %% (e.g) --threshold=0.1 => print orbital with more than 0.1 %% contribution")
    return parser.parse_args()


def is_this_row_for_coefficients(words: "list[str]") -> bool:
    if 8 <= len(words) <= 9:
        return True
    else:
        return False


def need_to_skip_this_line(words: "list[str]", is_reading_coefficients: bool) -> bool:
    if is_reading_coefficients and len(words) == 0:
        return False
    elif len(words) <= 1:
        return True
    else:
        return False


def need_to_write_results(words: "list[str]", is_reading_coefficients: bool) -> bool:
    if is_reading_coefficients and len(words) == 0:
        return True
    else:
        return False


def need_to_get_mo_sym_type(words: "list[str]", is_reading_coefficients: bool) -> bool:
    if not is_reading_coefficients and len(words) == 3 and words[0] == "Fermion" and words[1] == "ircop":
        return True
    return False


def need_to_start_reading_coefficients(words: "list[str]", is_reading_coefficients: bool) -> bool:
    if not is_reading_coefficients and len(words) == 6 and words[1] == "Electronic":
        return True
    return False


def get_dirac_filename(args: "argparse.Namespace") -> str:
    if not args.file:
        sys.exit("ERROR: DIRAC output file is not given. Please use -f option.")
    return args.file


def get_atom_num(num_str: "str | None", elem: str) -> int:
    if num_str is None:
        return 1
    elif num_str not in elem:
        sys.exit(f"ERROR: {num_str} is not in {elem}. Please check the molecule specification (-mol option).")
    return int(num_str)


def space_separated_parsing(line: str) -> "list[str]":
    words = re.split(" +", line.rstrip("\n"))
    return [word for word in words if word != ""]


def parse_molecule_input(args: "argparse.Namespace") -> Atoms:
    """
    Nested functions to parse molecule input
    """

    def parse_atom(input_list: "list[str]", line: str) -> str:
        letter = "".join(filter(str.isalpha, line)) or None
        if letter not in atoms.elements:
            sys.exit(f"ERROR: {line} is invalid.")
        if letter in atoms.atom_types:
            sys.exit(f"ERROR: {letter} is duplicated.\nYour input: {''.join(input_list)}")
        return letter

    def parse_atom_num(line: str) -> int:
        num_str = "".join(filter(str.isdigit, line)) or None
        return get_atom_num(num_str, line)

    def parse_atom_and_the_number_of_atoms(split_list: "list[str]") -> Atoms:
        for elem in split_list:

            letter = parse_atom(split_list, elem)
            atoms.atom_types.append(letter)

            num = parse_atom_num(elem)
            atoms.atom_nums.append(num)
        return atoms

    """
    Main function to parse molecule input
    """
    atoms = Atoms()
    if not args.molecule:
        sys.exit("ERROR: Molecule is not specified. Please use -mol option. (e.g. -mol Cu2O)")
    split_per_atom = re.findall("[A-Z][^A-Z]*", args.molecule)
    return parse_atom_and_the_number_of_atoms(split_per_atom)


def get_and_add_coefficient(words: "list[str]", atoms: Atoms, coefficients: Coefficients) -> None:
    """
    Nested functions to get and add coefficient
    """

    def get_orbital_type() -> "tuple[str, str]":
        atom_type = ""
        symmetry_type = ""
        if len(words) == 9:  # Ag pattern
            """
            (e.g)   words[2] = "Ag"
                    words[3] = "O"
                    symmetry_type = "Ag"
                    atom_type = "O"
            """
            symmetry_type = words[2]
            atom_type = words[3]
        elif len(words) == 8:
            """
            (e.g)   words[2] = "B3gO"
                    symmetry_type = "B3g"
                    atom_type = "O"
            """
            symmetry_type = words[2][:3]
            atom_type = words[2][3:]
        return symmetry_type, atom_type

    def add_orbital_type(atom_type: str, atom_orb_type: str) -> None:
        if atom_orb_type not in atoms.orbital_types:
            atoms.orbital_types.append(atom_orb_type)
            atoms.atom_list.append(atom_type)
            coefficients.mo_coefficient_list.append(0.0)

    def trim_words() -> "list[str]":
        if len(words) == 9:
            return words[4:]
        else:
            return words[3:]

    def get_coefficient() -> float:
        """
        (e.g)
        words = ["g400", "0.0000278056", "0.0000000000", "0.0000000000", "0.0000000000"]
        """
        alpha1: float = float(words[1])
        alpha2: float = float(words[2])
        beta1: float = float(words[3])
        beta2: float = float(words[4])
        return alpha1**2 + alpha2**2 + beta1**2 + beta2**2

    def check_atom_type() -> None:
        if atom_type not in atoms.atom_types:
            sys.exit(f"ERROR: atom type {atom_type} is not defined. Please check your --molecule or -mol option.")

    def add_coefficient(coefficient: float, atom_orb_type: str) -> None:
        magnification = atoms.atom_nums[atoms.atom_types.index(atom_type)]

        coefficient = magnification * coefficient
        coefficients.norm_const_sum += coefficient

        orb_type_idx = atoms.orbital_types.index(atom_orb_type)
        coefficients.mo_coefficient_list[orb_type_idx] += coefficient

    """
    Main function to get and add coefficient
    """
    symmetry_type, atom_type = get_orbital_type()
    # Trimming simplifies the structure of the function to get the coefficients.
    words = trim_words()
    coefficient = get_coefficient()

    orbital_type: str = words[0]  # orbital type (e.g. dxy)
    atom_orb_type = symmetry_type + atom_type + orbital_type

    add_orbital_type(atom_type, atom_orb_type)
    check_atom_type()
    add_coefficient(coefficient, atom_orb_type)

    return None


def write_results(
    atoms: Atoms,
    coefficients: Coefficients,
    threshold: float,
) -> None:
    """
    Nested functions to write results
    """

    def create_data_per_orbital_types(data_per_orb_types: Data_per_orbital_types) -> None:
        for (orb, coefficient, atom) in zip(atoms.orbital_types, coefficients.mo_coefficient_list, atoms.atom_list):
            atom_num = atoms.atom_nums[atoms.atom_types.index(atom)]
            data_per_orb_types.orbital_type.append(orb)
            data_per_orb_types.atom.append(atom)
            data_per_orb_types.mo_percentage.append(coefficient * 100 / (coefficients.norm_const_sum * atom_num))
        return

    def calculate_sum_of_mo_coefficient() -> float:
        return (sum([c for c in coefficients.mo_coefficient_list])) / coefficients.norm_const_sum

    """
    Main function to write results
    """
    data_per_orbital_types = Data_per_orbital_types()
    data_per_orbital_types.reset()  # Initialize the data_per_orbital_types
    create_data_per_orbital_types(data_per_orbital_types)
    coefficients.sum_of_mo_coefficient = calculate_sum_of_mo_coefficient()
    for idx, mo_percentage in enumerate(data_per_orbital_types.mo_percentage):
        if mo_percentage > threshold:
            for _ in range(atoms.atom_nums[atoms.atom_types.index(data_per_orbital_types.atom[idx])]):
                orb_type = str(data_per_orbital_types.orbital_type[idx]).ljust(11, " ")
                print(f"{orb_type}: {data_per_orbital_types.mo_percentage[idx]}  \t%")
    print(f"Normalization constant is {coefficients.norm_const_sum}")
    print(f"sum of coefficient {coefficients.sum_of_mo_coefficient}\n")

    return None


def main() -> None:

    is_reading_coefficients: bool = False
    electron_number: int = 0
    mo_sym_type: str = ""
    coefficients = Coefficients()

    args: "argparse.Namespace" = parse_args()
    threshold: float = args.threshold if args.threshold else 0.1
    dirac_file: str = get_dirac_filename(args)
    atoms: Atoms = parse_molecule_input(args)

    with open(dirac_file) as f:
        for line in f:

            words: "list[str]" = space_separated_parsing(line)

            if need_to_skip_this_line(words, is_reading_coefficients):
                continue

            elif need_to_get_mo_sym_type(words, is_reading_coefficients):
                mo_sym_type = words[2]

            elif need_to_start_reading_coefficients(words, is_reading_coefficients):
                """
                (e.g.)
                words = ["*", "Electronic", "eigenvalue", "no.", "22:", "-2.8417809384721"]
                """
                is_reading_coefficients = True
                electron_number = int(words[4][:-1])
                mo_energy = float(words[5])
                print("Electronic no.", electron_number, mo_sym_type, mo_energy)

            elif need_to_write_results(words, is_reading_coefficients):
                is_reading_coefficients = False
                write_results(atoms, coefficients, threshold)
                # Reset variables
                coefficients.reset()
                atoms.reset()

            elif is_reading_coefficients:
                if not is_this_row_for_coefficients(words):
                    continue
                get_and_add_coefficient(words, atoms, coefficients)


if __name__ == "__main__":
    main()
