import concurrent.futures
from enum import Enum, auto
from typing import Dict, List, Tuple

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.atoms import AtomInfo
from sum_dirac_dfcoef.coefficient import get_coefficient
from sum_dirac_dfcoef.data import DataAllMO, DataMO
from sum_dirac_dfcoef.eigenvalues import Eigenvalues
from sum_dirac_dfcoef.functions_info import FunctionsInfo
from sum_dirac_dfcoef.utils import debug_print, fast_deepcopy_pickle, space_separated_parsing


class STAGE(Enum):
    #                                                                             SKIP_READING_COEF
    #                                                                                    ↑↓
    # STAGE TRANSITION: INIT -> SKIP_AFTER_VECTOR_PRINT_LINE -> VECTOR_PRINT -> WAIT_END_READING_COEF -> END
    #                                                                             ↓               ↑
    #                                                                      WAIT_FIRST_COEF -> READING_COEF
    INIT = auto()
    SKIP_AFTER_VECTOR_PRINT_LINE = auto()
    VECTOR_PRINT = auto()
    WAIT_END_READING_COEF = auto()
    SKIP_READING_COEF = auto()
    WAIT_FIRST_COEF = auto()
    READING_COEF = auto()
    END = auto()


class PrivecProcessor:
    """This class has methods to read coefficients from the output file of DIRAC and store them in self.data_all_mo (final result).

    Attributes:
        dirac_output (List[str]): Output file strings of DIRAC
        stage (STAGE): Stage of reading coefficients
        is_electronic (bool): True if the current MO is electronic
        is_less_than_dirac_21 (bool): True if the current version of DIRAC is less than 21.
                                      Use this variable to determine which attribute (func_idx_dirac21/19) of AtomInfo should be used.
        eigenvalues (Eigenvalues): Eigenvalues
        mo_sym_type (str): Symmetry type of the current MO
        functions_info (FunctionsInfo): FunctionsInfo (used to handle which functions are used in the current MO)
        data_mo (DataMO): DataMO (temporary result of current reading MO)
        data_all_mo (DataAllMO): DataAllMO (final result)
        used_atom_info (Dict[str, AtomInfo]): Used AtomInfo
        current_atom_info (AtomInfo): Current AtomInfo
    """

    def __init__(self, dirac_output: List[str], functions_info: FunctionsInfo, eigenvalues: Eigenvalues) -> None:
        self.dirac_output = dirac_output
        self.stage = STAGE.INIT
        self.is_electronic = False
        self.is_less_than_dirac_21 = False
        self.eigenvalues = eigenvalues
        self.mo_sym_type = ""
        self.functions_info = functions_info
        self.data_mo = DataMO()
        self.data_all_mo = DataAllMO()
        self.used_atom_info: Dict[str, AtomInfo] = {}
        self.current_atom_info = AtomInfo()

    def read_privec_data(self, rank: int = 0) -> Tuple[DataAllMO, Eigenvalues]:
        """Read coefficients from the output file of DIRAC and store them in data_all_mo.

        self.data_all is the final result of this function. You can get all results from this variable except header information.
        """
        self.stage = STAGE.INIT
        mo_cnt = 0
        for line_str in self.dirac_output:
            if self.stage == STAGE.END:
                break  # End of reading coefficients
            elif self.stage == STAGE.SKIP_READING_COEF:
                if len(line_str.strip()) == 0:
                    self.transition_stage(STAGE.WAIT_END_READING_COEF)
                continue
            elif self.stage == STAGE.SKIP_AFTER_VECTOR_PRINT_LINE:
                self.transition_stage(STAGE.VECTOR_PRINT)
                continue

            words = space_separated_parsing(line_str)

            if self.need_to_skip_this_line(words):
                if self.stage == STAGE.VECTOR_PRINT and self.detect_next_titler(line_str):
                    msg = "WARNING: The next title is detected before the end of reading coefficients.\n\
In order to force DIRAC to print vector print, please add .ANALYZE and .PRIVEC option to the input file of DIRAC."
                    print(msg)
                    self.transition_stage(STAGE.END)
                elif self.stage == STAGE.READING_COEF:
                    if self.need_to_create_results_for_current_mo(words):
                        self.add_current_mo_data_to_data_all_mo()
                        self.transition_stage(STAGE.WAIT_END_READING_COEF)

            elif self.stage == STAGE.INIT:
                if self.check_start_vector_print(words):
                    self.transition_stage(STAGE.SKIP_AFTER_VECTOR_PRINT_LINE)

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
                    if args.parallel == 1 or mo_cnt % args.parallel == rank:
                        # Need to read coefficients of the current MO
                        self.add_coefficient(line_str)
                        self.transition_stage(STAGE.READING_COEF)
                    else:
                        # Don't need to read coefficients of the current MO because it is read by another process.
                        # multi-process version only
                        self.transition_stage(STAGE.SKIP_READING_COEF)
                    mo_cnt += 1

            elif self.stage == STAGE.READING_COEF:
                if self.is_this_row_for_coefficients(words):
                    self.add_coefficient(line_str)

        return self.data_all_mo, self.eigenvalues

    def transition_stage(self, new_stage: STAGE) -> None:
        self.stage = new_stage

    def is_this_row_for_coefficients(self, words: List[str]) -> bool:
        # min: 4 coefficients and other words => 5 words
        return True if 5 <= len(words) <= 9 and words[0].isdigit() else False

    def need_to_skip_this_line(self, words: List[str]) -> bool:
        return True if len(words) <= 1 else False

    def need_to_create_results_for_current_mo(self, words: List[str]) -> bool:
        return True if self.stage == STAGE.READING_COEF and len(words) <= 1 else False

    def detect_next_titler(self, line_str: str) -> bool:
        """Detect all characters are asterisk or space or line break."""
        stripped_line = line_str.strip().strip("\r\n")
        if len(stripped_line) == 0:
            return False
        return all(c == "*" for c in stripped_line)

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
        # https://github.com/RQC-HU/sum_dirac_dfcoef/issues/7#issuecomment-1377969626
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
        self.data_mo.sym_type = self.mo_sym_type
        self.data_mo.mo_info = mo_info
        self.used_atom_info.clear()  # reset used_atom_info because we need to delete used_atom_info of the previous MO
        self.current_atom_info = AtomInfo()  # reset current_atom_info because self.current_atom_info.count_remaining_functions() may be larger than 0

    def add_coefficient(self, line_str: str) -> None:
        def get_func_no() -> int:
            try:
                num_functions = int(line_str[:10])
            except ValueError as e:
                msg = f"num_functions must be integer, but got {line_str[:10]}"
                raise ValueError(msg) from e
            return num_functions

        def is_dirac_version_less_than_21():
            if self.is_less_than_dirac_21:
                return True
            for atom_info in self.functions_info[component_func][symmetry_label][atom_label].values():
                if atom_info.func_idx_dirac21.first <= num_functions <= atom_info.func_idx_dirac21.last:
                    # If the num_functions is include in the range of the current atom_info indices, it means that the current atom is DIRAC 21 or later.
                    return False
            return True

        def need_to_update_current_atom_info():
            first = self.current_atom_info.func_idx_dirac19.first if self.is_less_than_dirac_21 else self.current_atom_info.func_idx_dirac21.first
            last = self.current_atom_info.func_idx_dirac19.last if self.is_less_than_dirac_21 else self.current_atom_info.func_idx_dirac21.last
            if first <= num_functions <= last:
                return False
            return True

        num_functions = get_func_no()
        component_func = "large" if line_str[10] == "L" else ("small" if line_str[10] == "S" else "")  # CLS
        symmetry_label = line_str[12:15].strip()  # REP (e.g. "Ag "), symmetry_label="Ag"
        atom_label = line_str[15:18].strip()  # NAMN (e.g. "Cm "), atom_labe="Cm"
        gto_type = line_str[18:22].strip()  # GTOTYP (e.g. "s   "), gto_type="s"
        label = symmetry_label + atom_label

        if need_to_update_current_atom_info():
            self.is_less_than_dirac_21 = is_dirac_version_less_than_21()
            found_next_atom_info = False
            for atom_info in self.functions_info[component_func][symmetry_label][atom_label].values():
                start = atom_info.func_idx_dirac19.first if self.is_less_than_dirac_21 else atom_info.func_idx_dirac21.first
                end = atom_info.func_idx_dirac19.last if self.is_less_than_dirac_21 else atom_info.func_idx_dirac21.last
                if start <= num_functions <= end:
                    found_next_atom_info = True
                    self.current_atom_info = fast_deepcopy_pickle(atom_info)
                    self.used_atom_info[label] = atom_info
                    break
            if not found_next_atom_info:
                msg = f"The corresponding atom_info is not found in functions_info[{component_func}][{symmetry_label}][{atom_label}],\
                    list[atom_info] = {list(self.functions_info[component_func][symmetry_label][atom_label].values())}"
                raise Exception(msg)

        self.current_atom_info.decrement_function(gto_type)
        coef = get_coefficient(line_str, self.functions_info, self.current_atom_info.idx_within_same_atom)
        for idx in range(coef.multiplication):
            atom_idx = coef.idx_within_same_atom + idx
            self.data_mo.add_coefficient(atom_idx, fast_deepcopy_pickle(coef))

    def add_current_mo_data_to_data_all_mo(self) -> None:
        self.data_mo.fileter_coefficients_by_threshold()
        # add current MO data to data_all_mo
        # create a new DataMO object using pickle
        copy_data_mo = fast_deepcopy_pickle(self.data_mo)
        if self.is_electronic:
            self.data_all_mo.electronic.append(copy_data_mo)
            cur_sym = self.mo_sym_type
            if args.for_generator:
                self.eigenvalues.energies_used[cur_sym][self.data_mo.eigenvalue_no] = True
        else:
            self.data_all_mo.positronic.append(copy_data_mo)
        debug_print(f"End of reading {self.data_mo.eigenvalue_no}th MO")

    def read_privec_data_wrapper(self):
        """Read coefficients from the output file of DIRAC and store them in data_all_mo.
        This function is intended to wrap the processing of
        read_privec_data between the single-process and multiprocess versions.
        """

        def merge_energies_used(eigenvalues: Eigenvalues) -> None:
            """If multi-process version is used, eigenvalues.energies_used is stored in different memory space.
            Therefore, we need to merge eigenvalues.energies_used of each process.

            Args:
                eigenvalues (Eigenvalues): eigenvalues.energies_used of each process is stored in this variable.
            """
            for sym_type_key, val in eigenvalues.energies_used.items():
                for eigenvalue_no, is_found in val.items():
                    if is_found:
                        self.eigenvalues.energies_used[sym_type_key][eigenvalue_no] = True

        num_processes: int = args.parallel
        if num_processes > 1:
            # Multi-process version
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = [executor.submit(self.read_privec_data, i) for i in range(num_processes)]
                result_list = [future.result() for future in concurrent.futures.as_completed(futures)]
            result: Tuple[DataAllMO, Eigenvalues]
            for result in result_list:
                data_all_mo = result[0]
                self.data_all_mo.electronic.extend(data_all_mo.electronic)
                self.data_all_mo.positronic.extend(data_all_mo.positronic)
                eigenvalues = result[1]
                merge_energies_used(eigenvalues)
        else:
            # Single-process version
            self.read_privec_data()

        if args.for_generator:
            self.fill_non_moltra_range_electronic_eigenvalues()
        self.data_all_mo.sort_mo_sym_type()

    def fill_non_moltra_range_electronic_eigenvalues(self):
        self.is_electronic = True
        for sym_type_key, val in self.eigenvalues.energies_used.items():
            self.mo_sym_type = sym_type_key
            for eigenvalue_no, used in val.items():
                if not used:
                    self.data_mo.reset()
                    self.data_mo.eigenvalue_no = eigenvalue_no
                    self.data_mo.mo_info = self.get_mo_info(eigenvalue_no)
                    self.data_mo.sym_type = sym_type_key
                    self.data_mo.mo_energy = self.eigenvalues.energies[sym_type_key][eigenvalue_no]
                    copy_data_mo = fast_deepcopy_pickle(self.data_mo)
                    self.data_all_mo.electronic.append(copy_data_mo)
