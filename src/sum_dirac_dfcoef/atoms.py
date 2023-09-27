from typing import Annotated


from annotated_types import MaxLen
from pydantic import BaseModel


class AtomicOrbital(BaseModel, validate_assignment=True):
    atom: str = ""
    subshell: Annotated[str, MaxLen(max_length=1)] = "s"
    gto_type: str = "s"

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
    function_types: "set[str]" = set()

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
