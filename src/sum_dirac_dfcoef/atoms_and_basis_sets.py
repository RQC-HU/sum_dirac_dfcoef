from io import TextIOWrapper
import sys

from .utils import space_separated_parsing


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
