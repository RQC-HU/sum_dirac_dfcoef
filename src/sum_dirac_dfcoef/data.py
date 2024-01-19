from typing import ClassVar, Dict, List

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.coefficient import Coefficient


class DataMO:
    norm_const_sum: float = 0.0
    mo_energy: float = 0.0
    mo_info: str = ""
    sym_type: str = ""
    eigenvalue_no: int = 0
    coef_dict: ClassVar[Dict[str, Coefficient]] = {}
    coef_list: ClassVar[List[Coefficient]] = []

    def __repr__(self) -> str:
        return f"mo_info: {self.mo_info}, mo_energy: {self.mo_energy}, eigenvalue_no: {self.eigenvalue_no}, mo_sym_type: {self.sym_type}, coef_dict: {self.coef_dict}"

    def add_coefficient(self, coef: Coefficient) -> None:
        key = (coef.function_label, coef.start_idx)
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
    electronic: ClassVar[List[DataMO]] = []
    positronic: ClassVar[List[DataMO]] = []

    def __repr__(self) -> str:
        return f"electronic: {self.electronic}, positronic: {self.positronic}"

    def sort_mo_sym_type(self) -> None:
        self.electronic.sort(key=lambda mo: (mo.sym_type, mo.mo_energy))
        self.positronic.sort(key=lambda mo: (mo.sym_type, mo.mo_energy))

    def sort_mo_energy(self) -> None:
        self.electronic.sort(key=lambda mo: mo.mo_energy)
        self.positronic.sort(key=lambda mo: mo.mo_energy)
