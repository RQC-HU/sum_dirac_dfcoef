import copy
from enum import Enum, auto
from io import TextIOWrapper
from typing import List

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.atoms import AtomInfo
from sum_dirac_dfcoef.coefficient import get_coefficient
from sum_dirac_dfcoef.data import DataAllMO, DataMO
from sum_dirac_dfcoef.functions_info import FunctionsInfo
from sum_dirac_dfcoef.utils import (
    debug_print,
    space_separated_parsing,
)


class STAGE(Enum):
    # STAGE TRANSITION: INIT -> VECTOR_PRINT -> WAIT_END_READING_COEF ->
    # [MO_COEF -> READING_COEF -> WAIT_END_READING_COEF -> MO_COEF -> READING_COEF -> WAIT_END_READING_COEF -> ...] -> END
    INIT = auto()
    VECTOR_PRINT = auto()
    WAIT_END_READING_COEF = auto()
    MO_COEF = auto()
    READING_COEF = auto()
    END = auto()


class PrivecProcessor:
    def __init__(self, dirac_output: TextIOWrapper, functions_info: FunctionsInfo) -> None:
        self.dirac_output = dirac_output
        self.stage = STAGE.INIT
        self.is_electronic = False
        self.mo_sym_type = ""
        self.functions_info = functions_info
        self.data_mo = DataMO()
        self.data_all_mo = DataAllMO()
        self.used_atom_info: dict[str, AtomInfo] = {}
        self.current_atom_info = AtomInfo()

    def read_privec_data(self):
        """Read coefficients from the output file of DIRAC and store them in data_all_mo.

        self.data_all is the final result of this function. You can get all results from this variable except header information.
        """
        self.stage = STAGE.INIT
        for line_str in self.dirac_output:
            words = space_separated_parsing(line_str)

            if self.stage == STAGE.END:
                break  # End of reading coefficients

            elif self.need_to_skip_this_line(words):
                if self.stage == STAGE.READING_COEF:
                    if self.need_to_create_results_for_current_mo(words):
                        self.add_current_mo_data_to_data_all_mo()
                        self.transition_stage(STAGE.WAIT_END_READING_COEF)

            elif self.stage == STAGE.INIT:
                if self.check_start_vector_print(words):
                    self.transition_stage(STAGE.VECTOR_PRINT)

            elif self.stage == STAGE.VECTOR_PRINT:
                if self.need_to_get_mo_sym_type(words):
                    self.mo_sym_type = words[2]
                    self.transition_stage(STAGE.WAIT_END_READING_COEF)

            elif self.stage == STAGE.WAIT_END_READING_COEF:
                if self.need_to_get_mo_sym_type(words):
                    self.mo_sym_type = words[2]
                elif self.need_to_start_mo_section(words):
                    self.start_mo_section(words)
                    self.transition_stage(STAGE.MO_COEF)
                elif self.check_end_vector_print(words):
                    self.transition_stage(STAGE.END)

            elif self.stage == STAGE.MO_COEF:
                if self.is_this_row_for_coefficients(words):
                    self.add_coefficient(line_str)
                    self.transition_stage(STAGE.READING_COEF)

            elif self.stage == STAGE.READING_COEF:
                if self.is_this_row_for_coefficients(words):
                    self.add_coefficient(line_str)

        if not args.no_sort:
            self.data_all_mo.sort_mo_energy()

    def transition_stage(self, new_stage: STAGE) -> None:
        self.stage = new_stage

    def is_this_row_for_coefficients(self, words: List[str]) -> bool:
        # min: 4 coefficients and other words => 5 words
        return True if 5 <= len(words) <= 9 and words[0].isdigit() else False

    def need_to_skip_this_line(self, words: List[str]) -> bool:
        return True if len(words) <= 1 else False

    def need_to_create_results_for_current_mo(self, words: List[str]) -> bool:
        return True if self.stage == STAGE.READING_COEF and len(words) <= 1 else False

    def need_to_get_mo_sym_type(self, words: List[str]) -> bool:
        return True if len(words) == 3 and words[0] == "Fermion" and words[1] == "ircop" else False

    def need_to_start_mo_section(self, words: List[str]) -> bool:
        if words[1] in ("Electronic", "Positronic") and words[2] == "eigenvalue" and "no." in words[3]:
            return True
        return False

    def check_start_vector_print(self, words: List[str]) -> bool:
        # ****************************** Vector print ******************************
        if len(words) < 4:
            return False
        elif words[1] == "Vector" and words[2] == "print":
            return True
        return False

    def check_end_vector_print(self, words: List[str]) -> bool:
        # https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/issues/7#issuecomment-1377969626
        if len(words) >= 2 and self.stage == STAGE.WAIT_END_READING_COEF:
            return True
        return False

    def start_mo_section(self, words: List[str]) -> None:
        """
        (e.g.)
        words = ["*", "Electronic", "eigenvalue", "no.", "22:", "-2.8417809384721"]
        words = ["*", "Electronic", "eigenvalue", "no.122:", "-2.8417809384721"]
        """

        def set_is_electronic() -> None:
            if words[1] == "Positronic":
                self.is_electronic = False
            elif words[1] == "Electronic":
                self.is_electronic = True
            else:
                msg = f"ERROR: UnKnow MO type, MO_Type={words[1]}"
                raise ValueError(msg)

        def get_electron_num():
            try:
                electron_num = int(words[-2][:-1].replace("no.", ""))
            except ValueError:
                # If *** is printed, we have no information about what number this MO is.
                # Therefore, we assume that electron_num is the next number after prev_electron_num.
                prev_electron_num = self.data_mo.electron_num  # prev_electron is the number of electrons of the previous MO
                electron_num = prev_electron_num + 1
            return electron_num

        set_is_electronic()
        electron_num = get_electron_num()
        mo_energy = float(words[-1])
        mo_info = (
            f"{self.mo_sym_type} {electron_num}"
            if args.compress
            else (f"Electronic no. {electron_num} {self.mo_sym_type}" if self.is_electronic else f"Positronic no. {electron_num} {self.mo_sym_type}")
        )
        # Here is the start point of reading coefficients of the current MO
        self.data_mo.reset()  # reset data_mo because we need to delete data_mo of the previous MO
        self.data_mo.electron_num = electron_num
        self.data_mo.mo_energy = mo_energy
        self.data_mo.mo_info = mo_info
        self.used_atom_info.clear()  # reset used_atom_info because we need to delete used_atom_info of the previous MO

    def add_coefficient(self, line_str: str) -> None:
        component_func = "large" if line_str[10] == "L" else ("small" if line_str[10] == "S" else "")  # CLS
        symmetry_label = line_str[12:15].strip()  # REP (e.g. "Ag "), symmetry_label="Ag"
        atom_label = line_str[15:18].strip()  # NAMN (e.g. "Cm "), atom_labe="Cm"
        gto_type = line_str[18:22].strip()  # GTOTYP (e.g. "s   "), gto_type="s"
        label = symmetry_label + atom_label

        if self.current_atom_info.count_remaining_functions() == 0 or label != self.current_atom_info.label:
            # First, we need to read information about the current atom.
            if label not in self.used_atom_info:
                # It is the first time to read information about the current atom.
                cur_atom_start_idx = 1
            else:
                # It is not the first time to read information about the current atom.
                # So we need to read information about the previous atom from used_atom_info.
                # current start_idx = previous start_idx + previous mul
                cur_atom_start_idx = self.used_atom_info[label].start_idx + self.used_atom_info[label].mul
            # Validate start_idx
            if cur_atom_start_idx not in self.functions_info[component_func][symmetry_label][atom_label]:
                msg = f"start_idx={cur_atom_start_idx} is not found in functions_info[{component_func}][{symmetry_label}][{atom_label}]"
                raise Exception(msg)
            # We can get information about the current atom from functions_info with start_idx.
            self.current_atom_info = copy.deepcopy(self.functions_info[component_func][symmetry_label][atom_label][cur_atom_start_idx])
            # Update used_atom_info with current_atom_info
            self.used_atom_info[label] = copy.deepcopy(self.current_atom_info)

        self.current_atom_info.decrement_function(gto_type)
        self.data_mo.add_coefficient(get_coefficient(line_str, self.functions_info, self.current_atom_info.start_idx))

    def add_current_mo_data_to_data_all_mo(self) -> None:
        self.data_mo.fileter_coefficients_by_threshold()
        if self.is_electronic:
            self.data_all_mo.electronic.append(copy.deepcopy(self.data_mo))
        else:
            self.data_all_mo.positronic.append(copy.deepcopy(self.data_mo))
        debug_print(f"End of reading {self.data_mo.electron_num}th MO")
