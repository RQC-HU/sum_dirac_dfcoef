import copy
from enum import Enum, auto
from io import TextIOWrapper
from typing import Dict, List

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.atoms import AtomInfo
from sum_dirac_dfcoef.coefficient import get_coefficient
from sum_dirac_dfcoef.data import DataAllMO, DataMO
from sum_dirac_dfcoef.eigenvalues import Eigenvalues
from sum_dirac_dfcoef.functions_info import FunctionsInfo
from sum_dirac_dfcoef.utils import debug_print, space_separated_parsing


class STAGE(Enum):
    # STAGE TRANSITION: INIT -> VECTOR_PRINT -> WAIT_END_READING_COEF -> END
    #                                             ↓               ↑
    #                                      WAIT_FIRST_COEF -> READING_COEF
    INIT = auto()
    VECTOR_PRINT = auto()
    WAIT_END_READING_COEF = auto()
    WAIT_FIRST_COEF = auto()
    READING_COEF = auto()
    END = auto()


class PrivecProcessor:
    def __init__(self, dirac_output: TextIOWrapper, functions_info: FunctionsInfo, eigenvalues: Eigenvalues) -> None:
        self.dirac_output = dirac_output
        self.stage = STAGE.INIT
        self.is_electronic = False
        self.eigenvalues = eigenvalues
        self.mo_sym_type = ""
        self.functions_info = functions_info
        self.data_mo = DataMO()
        self.data_all_mo = DataAllMO()
        self.used_atom_info: Dict[str, AtomInfo] = {}
        self.current_atom_info = AtomInfo()

    def read_privec_data(self):
        """Read coefficients from the output file of DIRAC and store them in data_all_mo.

        self.data_all is the final result of this function. You can get all results from this variable except header information.
        """
        self.stage = STAGE.INIT
        for line_str in self.dirac_output:
            words = space_separated_parsing(line_str)

            if self.stage == STAGE.END:
                if args.for_generator:
                    self.fill_non_moltra_range_electronic_eigenvalues()
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
                    self.transition_stage(STAGE.WAIT_FIRST_COEF)
                elif self.check_end_vector_print(words):
                    self.transition_stage(STAGE.END)

            elif self.stage == STAGE.WAIT_FIRST_COEF:
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

    def get_mo_info(self, eigenvalue_no: int) -> str:
        if args.compress:
            return f"{self.mo_sym_type} {eigenvalue_no}"
        elif self.is_electronic:
            return f"Electronic no. {eigenvalue_no} {self.mo_sym_type}"
        else:
            return f"Positronic no. {eigenvalue_no} {self.mo_sym_type}"

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

        def get_eigenvalue_no():
            try:
                eigenvalue_no = int(words[-2][:-1].replace("no.", ""))
            except ValueError:
                # If *** is printed, we have no information about what number this MO is.
                # Therefore, we assume that eigenvalue_no is the next number after prev_eigenvalue_no.
                prev_eigenvalue_no = self.data_mo.eigenvalue_no  # prev_electron is the number of electrons of the previous MO
                eigenvalue_no = prev_eigenvalue_no + 1
            return eigenvalue_no

        set_is_electronic()
        eigenvalue_no = get_eigenvalue_no()
        mo_energy = float(words[-1])
        mo_info = self.get_mo_info(eigenvalue_no)

        # Here is the start point of reading coefficients of the current MO
        self.data_mo.reset()  # reset data_mo because we need to delete data_mo of the previous MO
        self.data_mo.eigenvalue_no = eigenvalue_no
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
            cur_sym = self.mo_sym_type
            if args.for_generator:
                self.eigenvalues.energies_used[cur_sym][self.data_mo.eigenvalue_no] = True
        else:
            self.data_all_mo.positronic.append(copy.deepcopy(self.data_mo))
        debug_print(f"End of reading {self.data_mo.eigenvalue_no}th MO")

    def fill_non_moltra_range_electronic_eigenvalues(self):
        self.is_electronic = True
        print("Filling non-moltra range electronic eigenvalues...")
        print(self.eigenvalues.energies_used)
        for sym_type_key, val in self.eigenvalues.energies_used.items():
            self.mo_sym_type = sym_type_key
            for eigenvalue_no, is_found in val.items():
                if not is_found:
                    self.data_mo.reset()
                    self.data_mo.eigenvalue_no = eigenvalue_no
                    self.data_mo.mo_info = self.get_mo_info(eigenvalue_no)
                    self.data_mo.mo_energy = self.eigenvalues.energies[sym_type_key][eigenvalue_no - 1]
                    self.data_all_mo.electronic.append(copy.deepcopy(self.data_mo))
