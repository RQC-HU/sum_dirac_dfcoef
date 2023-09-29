from collections import OrderedDict
from typing import Set, OrderedDict as ODict


from pydantic import BaseModel, validator


class AtomInfo:
    start_idx: int
    label: str
    mul: int
    functions: ODict[str, int]

    def __init__(self, start_idx: int = 0, label: str = "", multiplicity: int = 0) -> None:
        self.start_idx = start_idx
        self.label = label
        self.mul = multiplicity
        self.functions = OrderedDict()

    def __repr__(self) -> str:
        return f"start_idx: {self.start_idx}, mul: {self.mul}, functions: {self.functions}"

    def add_function(self, gto_type: str, num_functions: int) -> None:
        self.functions[gto_type] = num_functions

    def decrement_function(self, gto_type: str) -> None:
        try:
            self.functions[gto_type] -= 1
            if self.functions[gto_type] < 0:
                raise ValueError(f"Too many functions detected. self.functions[{gto_type}] must be >= 0, but got {self.functions[gto_type]}")
        except KeyError:
            raise KeyError(f"{gto_type} is not in current_atom_info: {self.functions.keys()}")

    def count_remaining_functions(self) -> int:
        return sum(self.functions.values())

    def get_remaining_functions(self) -> "ODict[str, int]":
        return OrderedDict({k: v for k, v in self.functions.items() if v > 0})


class AtomicOrbital(BaseModel, validate_assignment=True):
    atom: str = ""
    subshell: str = "s"
    gto_type: str = "s"

    @validator("subshell")
    def validate_subshell(cls, v: str) -> str:
        order_of_subshell = "spdfghiklmnoqrtuvwxyz"
        if v not in order_of_subshell:
            raise ValueError(f"subshell must be one of '{order_of_subshell}', but got '{v}'")
        if len(v) != 1:
            raise ValueError("subshell must be one character")
        return v

    def reset(self):
        self.atom = ""
        self.subshell = "s"
        self.gto_type = "s"

    def set(self, atom: str, subshell: str, gto_type: str):
        self.atom = atom
        self.subshell = subshell
        self.gto_type = gto_type


class AtomicOrbitals(BaseModel, validate_assignment=True):
    prev_ao: AtomicOrbital = AtomicOrbital()
    current_ao: AtomicOrbital = AtomicOrbital()
    start_idx: int = 1
    function_types: Set[str] = set()

    def reset(self):
        self.prev_ao.reset()
        self.current_ao.reset()
        self.start_idx = 1
        self.function_types.clear()


def is_different_atom(ao: AtomicOrbitals, function_label: str) -> bool:
    def is_reverse_subshell() -> bool:
        order_of_subshell = "spdfghiklmnoqrtuvwxyz"
        if order_of_subshell.index(ao.prev_ao.subshell) > order_of_subshell.index(ao.current_ao.subshell):
            return True
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
