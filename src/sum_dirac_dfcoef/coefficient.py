from typing import Dict

from pydantic import BaseModel


class Coefficient(BaseModel, validate_assignment=True):
    vector_num: int
    function_label: str
    need_identifier: bool
    coefficient: float
    start_idx: int
    multiplication: int

    def __repr__(self) -> str:
        super().__repr__()
        return f"vector_num: {self.vector_num}, function_label: {self.function_label}, coefficient: {self.coefficient}, start_idx: {self.start_idx}, multiplication: {self.multiplication}"


class Coefficients:
    norm_const_sum: float = 0.0
    coef_dict: Dict[str, Coefficient] = dict()
    coef_list: "list[Coefficient]" = list()
    mo_energy: float = 0.0
    mo_info: str = ""

    def __repr__(self) -> str:
        # return f"norm_const_sum: {self.norm_const_sum}, coef_list: {self.coef_list}"
        return f"norm_const_sum: {self.norm_const_sum}, coef: {[coef.coefficient for coef in self.coef_list]}"

    def add_coefficient(self, coef: Coefficient) -> None:
        # self.coef_list.append(coef)
        if coef.function_label in self.coef_dict:
            self.coef_dict[coef.function_label].coefficient += coef.coefficient
        else:
            self.coef_dict[coef.function_label] = coef
        self.norm_const_sum += coef.coefficient * coef.multiplication

    def reset(self):
        self.norm_const_sum = 0.0
        self.mo_energy = 0.0
        self.mo_info = ""
        self.coef_dict.clear()
        self.coef_list.clear()
