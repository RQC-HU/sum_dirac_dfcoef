from collections import OrderedDict
from typing import OrderedDict as ODict
from typing import Set

from pydantic import BaseModel, validator

from sum_dirac_dfcoef.subshell import subshell_order


class LastFuncNum:
    gerade: int = 0  # Number of functions with gerade symmetry
    ungerade: int = 0  # Number of functions with ungerade symmetry
    no_inv_sym: int = 0  # Number of functions with no inversion symmetry
    sum_all: int = 0  # Sum of all functions

    def add(self, symmetry: str, num_functions: int) -> None:
        if symmetry[-1].lower() == "g":
            self.gerade += num_functions
        elif symmetry[-1].lower() == "u":
            self.ungerade += num_functions
        else:
            self.no_inv_sym += num_functions
        self.sum_all += num_functions

    def get(self, symmetry: str) -> int:
        if symmetry[-1].lower() == "g":
            return self.gerade
        elif symmetry[-1].lower() == "u":
            return self.ungerade
        else:
            return self.no_inv_sym

    def get_all(self) -> int:
        return self.sum_all


class FuncIndices:
    first: int
    last: int

    def __init__(self, first: int = 0, last: int = 0) -> None:
        self.first = first
        self.last = last


class AtomInfo:
    idx_within_same_atom: int
    label: str
    mul: int
    functions: ODict[str, int]
    func_idx_per_sym: FuncIndices  # For DIRAC < 21
    func_idx_all_sym: FuncIndices  # For DIRAC >= 21
    last_func_num: int

    def __init__(
        self,
        idx_within_same_atom: int = 0,
        label: str = "",
        multiplicity: int = 0,
        last_func_num: int = 0,
        func_idx_per_sym: FuncIndices = None,
        func_idx_all_sym: FuncIndices = None,
    ) -> None:
        self.idx_within_same_atom = idx_within_same_atom
        self.label = label
        self.mul = multiplicity
        self.functions = OrderedDict()
        self.func_idx_per_sym = func_idx_per_sym if func_idx_per_sym is not None else FuncIndices()
        self.func_idx_all_sym = func_idx_all_sym if func_idx_all_sym is not None else FuncIndices()
        self.last_func_num = last_func_num

    def __repr__(self) -> str:
        return f"idx_within_same_atom: {self.idx_within_same_atom}, mul: {self.mul}, last_func_num: {self.last_func_num}, func_idx_per_sym.first: {self.func_idx_per_sym.first}, func_idx_per_sym.last: {self.func_idx_per_sym.last}, func_idx_all_sym.first: {self.func_idx_all_sym.first}, func_idx_all_sym.last: {self.func_idx_all_sym.last}, functions: {self.functions}"

    def add_function(self, gto_type: str, num_functions: int) -> None:
        self.functions[gto_type] = num_functions

    def update_last_func_num(self, symmetry: str, last_func_num: LastFuncNum) -> None:
        self.last_func_num = last_func_num.get(symmetry)
        self.func_idx_per_sym.last = self.last_func_num
        # self.func_idx_all_sym.last = last_func_num.get_all()

    def decrement_function(self, gto_type: str) -> None:
        try:
            self.functions[gto_type] -= 1
            if self.functions[gto_type] < 0:
                msg = f"Too many functions detected. self.functions[{gto_type}] must be >= 0, but got {self.functions[gto_type]}"
                raise ValueError(msg)
        except KeyError as e:
            msg = f"self.functions[{gto_type}] is not found in self.functions: {self.functions.keys()}"
            raise KeyError(msg) from e


    def get_remaining_functions(self) -> "ODict[str, int]":
        return OrderedDict({k: v for k, v in self.functions.items() if v > 0})


class AtomicOrbital(BaseModel, validate_assignment=True):
    atom: str = ""
    subshell: str = "s"
    gto_type: str = "s"

    @validator("subshell")
    def validate_subshell(cls, v: str) -> str:  # noqa: N805 (pydantic method)
        if v not in subshell_order.subshell_order:
            msg = f"subshell must be one of '{subshell_order.subshell_order}', but got '{v}'"
            raise ValueError(msg)
        if len(v) != 1:
            msg = f"subshell must be one character, but got '{v}'"
            raise ValueError(msg)
        return v

    def reset(self):
        self.atom = ""
        self.subshell = "s"
        self.gto_type = "s"

    def update(self, atom: str, subshell: str, gto_type: str):
        self.atom = atom
        self.subshell = subshell
        self.gto_type = gto_type


class AtomicOrbitals(BaseModel, validate_assignment=True):
    prev_ao: AtomicOrbital = AtomicOrbital()
    current_ao: AtomicOrbital = AtomicOrbital()
    idx_within_same_atom: int = 1
    function_types: Set[str] = set()

    def reset(self):
        self.prev_ao.reset()
        self.current_ao.reset()
        self.idx_within_same_atom = 1
        self.function_types.clear()


def is_different_atom(ao: AtomicOrbitals, function_label: str) -> bool:
    def is_reverse_subshell() -> bool:
        prev_subshell_idx = subshell_order.subshell_order.index(ao.prev_ao.subshell)
        current_subshell_idx = subshell_order.subshell_order.index(ao.current_ao.subshell)
        if prev_subshell_idx > current_subshell_idx:
            return True
        elif prev_subshell_idx < current_subshell_idx:
            return False
        else:  # Same subshell
            subshell_idx = subshell_order.subshell_order.index(ao.prev_ao.subshell)
            prev_gto_idx = subshell_order.gto_label_order[subshell_idx].index(ao.prev_ao.gto_type)
            current_gto_idx = subshell_order.gto_label_order[subshell_idx].index(ao.current_ao.gto_type)
            if prev_gto_idx > current_gto_idx:  # reverse subshell. e.g. ao.prev_ato.gto_type = "pz", ao.current_ao.gto_type = "px"
                return True
            else:
                return False

    if ao.prev_ao.atom != ao.current_ao.atom:
        return True
    elif function_label in ao.function_types:
        # They have the same element label but different atoms.
        # Because the function_label of an atom is combined in one line,
        # it is a different atom if the same function_label appears.
        return True
    elif is_reverse_subshell():
        # e.g. "C  p" -> "C  s"
        # Different atom
        return True
    return False  # Same atom


ao = AtomicOrbitals()
