from collections import OrderedDict
from typing import List, Optional
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.coefficient import Coefficient


class CoefKey:
    """This class is used to store the key of the coefficient dictionary.

    Attributes:
        atom_label (str): The atom label to identify the vector.
        azimuthal_label (str): The azimuthal quantum number label to identify the vector.
        need_identifier (bool): Whether print the atom_index or not.
                                This field is only for passing this flag to the file_writer::write_mo_data method.
        atom_index (int): The atom index to identify the vector.
        symmetry_label (str): The symmetry label to identify the vector.
        magnetic_label (str): The magnetic quantum number label to identify the vector.
    """

    atom_label: str
    azimuthal_label: str
    need_identifier: bool
    atom_idx: int
    symmetry_label: str
    magnetic_label: str

    def __init__(self, atom_idx: int, coef: Coefficient) -> None:
        self.atom_label = coef.atom_label
        self.azimuthal_label = coef.azimuthal_label
        self.need_identifier = coef.need_identifier
        self.atom_idx = -1 if args.ignore_atom_num else atom_idx
        self.symmetry_label = "" if args.ignore_sym else coef.symmetry_label
        self.magnetic_label = "" if args.ignore_ml else coef.magnetic_label

    def __repr__(self) -> str:
        return f"CoefKey(atom_label: {self.atom_label}, \
azimuthal_label: {self.azimuthal_label}, \
need_identifier: {self.need_identifier}, \
atom_idx: {self.atom_idx}, \
symmetry_label: {self.symmetry_label}, \
magnetic_label: {self.magnetic_label})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CoefKey):
            return NotImplemented
        return (
            self.atom_label == other.atom_label
            and self.azimuthal_label == other.azimuthal_label
            and self.atom_idx == other.atom_idx
            and self.symmetry_label == other.symmetry_label
            and self.magnetic_label == other.magnetic_label
        )

    def __hash__(self) -> int:
        return hash((self.atom_label, self.azimuthal_label, self.atom_idx, self.symmetry_label, self.magnetic_label))


class DataMO:
    """This class is used to store the specific MO coefficient information.

    Attributes:
        norm_const_sum (float): The sum of the coefficients of the MO.
        mo_energy (float): The energy of the MO.
        mo_info (str): The information of the MO. (e.g. "Electronic no. 1 E1g")
        sym_type (str): The symmetry type of the MO. (e.g. "E1g")
        eigenvalue_no (int): The eigenvalue number of the MO.
        coef_dict (ODict[CoefKey, float]): The dictionary of the coefficients of the MO.
                                           It is sorted before printing to the file by calling the filter_coefficients_by_threshold method.
    """

    norm_const_sum: float = 0.0
    mo_energy: float = 0.0
    mo_info: str = ""
    sym_type: str = ""
    eigenvalue_no: int = 0
    coef_dict: ODict[CoefKey, float]

    def __init__(
        self,
        mo_info: str = "",
        mo_energy: float = 0.0,
        eigenvalue_no: int = 0,
        sym_type: str = "",
        coef_dict: Optional[ODict[CoefKey, float]] = None,
    ) -> None:
        self.mo_info = mo_info
        self.mo_energy = mo_energy
        self.eigenvalue_no = eigenvalue_no
        self.sym_type = sym_type
        self.coef_dict = coef_dict if coef_dict is not None else OrderedDict()

    def __repr__(self) -> str:
        return f"DataMO(mo_info: {self.mo_info}, mo_energy: {self.mo_energy}, eigenvalue_no: {self.eigenvalue_no}, mo_sym_type: {self.sym_type}, coef_dict: {self.coef_dict})"

    def add_coefficient(self, atom_idx: int, coef: Coefficient) -> None:
        key = CoefKey(atom_idx, coef)
        if key in self.coef_dict:
            self.coef_dict[key] += coef.coefficient
        else:
            self.coef_dict[key] = coef.coefficient
        self.norm_const_sum += coef.coefficient

    def reset(self):
        self.norm_const_sum = 0.0
        self.mo_energy = 0.0
        self.mo_info = ""
        self.eigenvalue_no = 0
        self.coef_dict.clear()

    def filter_coefficients_by_threshold(self) -> None:
        self.coef_dict = OrderedDict((key, coef) for key, coef in self.coef_dict.items() if abs(coef / self.norm_const_sum * 100) >= args.threshold)
        self.coef_dict = OrderedDict(sorted(self.coef_dict.items(), key=lambda x: x[1], reverse=True))


class DataAllMO:
    """This class stores all electronic and positronic MO data.

    This class holds information used for output excluding header information.

    Attributes:
        electronic (List[DataMO]): The list of the electronic MO data.
        positronic (List[DataMO]): The list of the positronic MO data.
    """

    electronic: List[DataMO]
    positronic: List[DataMO]

    def __init__(self, electronic: Optional[List[DataMO]] = None, positronic: Optional[List[DataMO]] = None) -> None:
        self.electronic = electronic if electronic is not None else []
        self.positronic = positronic if positronic is not None else []

    def __repr__(self) -> str:
        return f"DataAllMO(electronic: {self.electronic}, positronic: {self.positronic})"

    def sort_mo_sym_type(self) -> None:
        self.electronic.sort(key=lambda mo: (mo.sym_type, mo.mo_energy))
        self.positronic.sort(key=lambda mo: (mo.sym_type, mo.mo_energy))

    def sort_mo_energy(self) -> None:
        self.electronic.sort(key=lambda mo: mo.mo_energy)
        self.positronic.sort(key=lambda mo: mo.mo_energy)
