from typing import Dict

from .args import args
from .coefficient import Coefficient


class Data_MO:
    norm_const_sum: float = 0.0
    coef_dict: Dict[str, Coefficient] = dict()
    coef_list: "list[Coefficient]" = list()
    electron_num: int = 0
    mo_energy: float = 0.0
    mo_info: str = ""

    def __repr__(self) -> str:
        # return f"norm_const_sum: {self.norm_const_sum}, coef_list: {self.coef_list}"
        return f"norm_const_sum: {self.norm_const_sum}, coef: {[coef.coefficient for coef in self.coef_list]}"

    def add_coefficient(self, coef: Coefficient) -> None:
        if coef.function_label in self.coef_dict:
            self.coef_dict[coef.function_label].coefficient += coef.coefficient
        else:
            self.coef_dict[coef.function_label] = coef
        self.norm_const_sum += coef.coefficient * coef.multiplication

    def reset(self):
        self.norm_const_sum = 0.0
        self.mo_energy = 0.0
        self.mo_info = ""
        self.electron_num = 0
        self.coef_dict.clear()
        self.coef_list.clear()

    def fileter_coefficients_by_threshold(self) -> None:
        self.coef_list = [coef for coef in self.coef_dict.values() if abs(coef.coefficient / self.norm_const_sum * 100) >= args.threshold]
        self.coef_list.sort(key=lambda coef: coef.coefficient, reverse=True)


class Data_All_MO:
    electronic: "list[Data_MO]" = []
    positronic: "list[Data_MO]" = []

    def __repr__(self) -> str:
        return f"electronic: {self.electronic}, positronic: {self.positronic}"
