from typing import Dict, List, Optional

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.coefficient import Coefficient


class DataMO:
    """ This class is used to store the specific MO coefficient information.

    Attributes:
        norm_const_sum (float): The sum of the coefficients of the MO.
        mo_energy (float): The energy of the MO.
        mo_info (str): The information of the MO. (e.g. "Electronic no. 1 E1g")
        sym_type (str): The symmetry type of the MO. (e.g. "E1g")
        eigenvalue_no (int): The eigenvalue number of the MO.
        coef_dict (Dict[str, Coefficient]): The dictionary of the coefficients, in the order written in the DIRAC output file.
        coef_list (List[Coefficient]): The list of the coefficients, which is sorted by the coefficient value.
    """
    norm_const_sum: float = 0.0
    mo_energy: float = 0.0
    mo_info: str = ""
    sym_type: str = ""
    eigenvalue_no: int = 0
    coef_dict: Dict[str, Coefficient]
    coef_list: List[Coefficient]

    def __init__(
        self,
        mo_info: str = "",
        mo_energy: float = 0.0,
        eigenvalue_no: int = 0,
        sym_type: str = "",
        coef_dict: Optional[Dict[str, Coefficient]] = None,
        coef_list: Optional[List[Coefficient]] = None,
    ) -> None:
        self.mo_info = mo_info
        self.mo_energy = mo_energy
        self.eigenvalue_no = eigenvalue_no
        self.sym_type = sym_type
        self.coef_dict = coef_dict if coef_dict is not None else {}
        self.coef_list = coef_list if coef_list is not None else []

    def __repr__(self) -> str:
        return f"DataMO(mo_info: {self.mo_info}, mo_energy: {self.mo_energy}, eigenvalue_no: {self.eigenvalue_no}, mo_sym_type: {self.sym_type}, coef_dict: {self.coef_dict})"

    def add_coefficient(self, coef: Coefficient) -> None:
        key = (coef.function_label, coef.idx_within_same_atom)
        if not args.ignore_ml:
            key += (coef.magnetic_label,)
        if key in self.coef_dict:
            self.coef_dict[key].coefficient += coef.coefficient
        else:
            self.coef_dict[key] = coef
        self.norm_const_sum += coef.coefficient * coef.multiplication

    def reset(self):
        self.norm_const_sum = 0.0
        self.mo_energy = 0.0
        self.mo_info = ""
        self.eigenvalue_no = 0
        self.coef_dict.clear()
        self.coef_list.clear()

    def fileter_coefficients_by_threshold(self) -> None:
        self.coef_list = [coef for coef in self.coef_dict.values() if abs(coef.coefficient / self.norm_const_sum * 100) >= args.threshold]
        self.coef_list.sort(key=lambda coef: coef.coefficient, reverse=True)


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
