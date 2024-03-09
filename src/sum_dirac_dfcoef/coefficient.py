from pydantic import BaseModel

from sum_dirac_dfcoef.functions_info import FunctionsInfo
from sum_dirac_dfcoef.utils import is_float, space_separated_parsing


class Coefficient(BaseModel, validate_assignment=True):
    """ This class is used to store the specific coefficient information.

    e.g.   1  L Ag U  s      -0.0012445973  0.0000000000  0.0000000000  0.0000000000
    -> Coefficient(vector_num=1, function_label="AgUs", need_identifier=False, coefficient=1.54902243916729E-6, idx_within_same_atom=1, multiplication=1)

    Attributes:
        vector_num (int): Serial number of the vector
        function_label (str): Information about the vector to identify the vector
        need_identifier (bool): Whether the idx_within_same_atom is needed to print to the output file or not.
        coefficient (float): Coefficient value
        idx_within_same_atom (int): Index of the order of atoms in the same AtomLabel
        multiplication (int): The multiplicity of the specified atom label.
    """
    vector_num: int
    function_label: str
    need_identifier: bool
    coefficient: float
    idx_within_same_atom: int
    multiplication: int

    def __repr__(self) -> str:
        super().__repr__()
        return f"Coefficient(vector_num: {self.vector_num}, \
function_label: {self.function_label}, \
need_identifier: {self.need_identifier}, \
coefficient: {self.coefficient}, \
idx_within_same_atom: {self.idx_within_same_atom}, \
multiplication: {self.multiplication})"


def get_coefficient(line_str: str, orbitals: FunctionsInfo, idx_within_same_atom: int) -> Coefficient:
    """
    This function parses the line that contains the coefficient and returns the Coefficient object.
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

    def parse_line() -> Coefficient:
        """
        This function parses the line that contains the coefficient.

        line write source: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L440-442
        line format source: https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirout.F#L453
        line format : FORMAT(3X,I5,2X,A12,2X,4F14.10)
        """
        words = space_separated_parsing(line_str)
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
        component_func = "large" if line_str[10] == "L" else ("small" if line_str[10] == "S" else "")  # CLS
        symmetry_label = line_str[12:15].strip()  # REP (e.g. "Ag ")
        atom_label = line_str[15:18].strip()  # NAMN (e.g. "Cm "), atom_labe="Cm"
        function_label = line_str[12:22].strip().replace(" ", "")  # REP + NAMN + GTOTYP (e.g. "Ag Cm s   " => "AgCms")

        # COEF (4F14.10)
        # coefficients = [line_str[24:38], line_str[38:52], line_str[52:66], line_str[66:80]]
        coef_num = 4
        coef_len = 14
        coef_start_idx = 24
        coefficient = sum(
            pow(float(line_str[i : i + coef_len]) if is_float(line_str[i : i + coef_len]) else -100, 2)
            for i in range(coef_start_idx, coef_start_idx + coef_len * coef_num, coef_len)
        )

        need_identifier = (
            True if len(orbitals[component_func][symmetry_label][atom_label]) > 1 or orbitals[component_func][symmetry_label][atom_label][idx_within_same_atom].mul > 1 else False
        )
        multiplication = int(orbitals[component_func][symmetry_label][atom_label][idx_within_same_atom].mul)

        return Coefficient(
            vector_num=vec_num,
            function_label=function_label,
            need_identifier=need_identifier,
            coefficient=coefficient,
            idx_within_same_atom=idx_within_same_atom,
            multiplication=multiplication,
        )

    """
    Main function to get coefficient
    """
    coef = parse_line()

    return coef
