import re
from bisect import bisect_left, bisect_right
from collections import OrderedDict
from io import TextIOWrapper
from typing import List
from typing import OrderedDict as ODict

from sum_dirac_dfcoef.data import DataAllMO, DataMO
from sum_dirac_dfcoef.eigenvalues import Eigenvalues
from sum_dirac_dfcoef.electron_num import get_electron_num_from_input, get_electron_num_from_scf_field
from sum_dirac_dfcoef.moltra import MoltraInfo


class HeaderInfo:
    """Class to store header information for the sum_dirac_dfcoef module.

    Attributes:
        header (str): Header for the sum_dirac_dfcoef module.
        subheader1 (str): First subheader for the sum_dirac_dfcoef module.
        subheader2 (str): Second subheader for the sum_dirac_dfcoef module.
    """

    def __init__(self):
        self.moltra_info = MoltraInfo()
        self.eigenvalues = Eigenvalues()
        self.electrons = 0

    def read_header_info(self, dirac_output: TextIOWrapper) -> None:
        """Read the header information from the output file of DIRAC

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
            None: class attributes are updated
        """
        dirac_output.seek(0)
        self.__read_electron_number(dirac_output)
        dirac_output.seek(0)
        self.__validate_eigpri_option(dirac_output)
        dirac_output.seek(0)
        self.__read_moltra(dirac_output)
        self.__read_eigenvalues(dirac_output)
        self.__duplicate_moltra_str()

    def __read_electron_number(self, dirac_output: TextIOWrapper) -> None:
        self.electrons = get_electron_num_from_input(dirac_output)
        if self.electrons == 0:
            self.electrons = get_electron_num_from_scf_field(dirac_output)

    def __validate_eigpri_option(self, dirac_output: TextIOWrapper) -> None:
        self.eigenvalues.validate_eigpri_option(dirac_output)

    def __read_eigenvalues(self, dirac_output: TextIOWrapper) -> None:
        self.eigenvalues.get_eigenvalues(dirac_output)

    def __read_moltra(self, dirac_output: TextIOWrapper) -> None:
        self.moltra_info.read_moltra_section(dirac_output)

    def __duplicate_moltra_str(self) -> None:
        # Duplicate the moltra range string if it is not enough
        if self.moltra_info.is_default:
            # Set the default range string
            for _ in range(len(self.eigenvalues.shell_num)):
                self.moltra_info.range_str.append("ENERGY -20.0 10.0 1.0")
        if len(self.moltra_info.range_str) != len(self.eigenvalues.shell_num):
            # one-line input is allowed, duplicate the first line
            # https://gitlab.com/dirac/dirac/-/blob/ea717cdb294035d8af3ebe2b1e00cf94f1c1a6b7/src/moltra/trainp.F#L592-600
            for _ in range(len(self.eigenvalues.shell_num) - len(self.moltra_info.range_str)):
                self.moltra_info.range_str.append(self.moltra_info.range_str[0])

    def calculate_moltra_idx_range(self, data_all_mo: DataAllMO) -> None:
        keys = list(self.eigenvalues.shell_num.keys())
        for i, item in enumerate(self.moltra_info.range_str):
            symmetry_type = keys[i]
            if "ALL" == item.upper():
                self.moltra_info.range_dict[symmetry_type] = f"1..{len(self.eigenvalues.energies[symmetry_type])}"
            elif "ENERGY" in item.upper():
                self.moltra_info.range_dict[symmetry_type] = self.__parse_energy_str(item, symmetry_type, data_all_mo)
            else:
                self.moltra_info.range_dict[symmetry_type] = self.__parse_range_str(item, symmetry_type)

    def __parse_energy_str(self, energy_str: str, symmetry_type: str, data_all_mo: DataAllMO) -> str:
        """Parse the energy string

        Args:
            energy_str (str): Energy string

        Returns:
            str: Range string
        """

        def get_min_energy_idx(min_energy: float, step: float) -> int:
            energies = self.eigenvalues.energies[symmetry_type]
            # Find the index of the minimum energy without exceeding the step
            cur_idx = bisect_left(energies, min_energy)
            cur_min_energy = energies[cur_idx]
            # Search for the minimum energy index
            while cur_idx > 0:
                next_energy = energies[cur_idx - 1]
                if abs(cur_min_energy - next_energy) > step:
                    break  # Found the minimum energy index
                cur_min_energy = next_energy
                cur_idx -= 1
            return cur_idx

        def get_max_energy_idx(max_energy: float, step: float) -> int:
            energies = self.eigenvalues.energies[symmetry_type]
            # Find the index of the maximum energy without exceeding the step
            cur_idx = bisect_right(energies, max_energy)
            cur_max_energy = energies[cur_idx - 1]
            while cur_idx < len(energies):
                next_energy = energies[cur_idx]
                if abs(next_energy - cur_max_energy) > step:
                    break  # Found the maximum energy index
                cur_max_energy = next_energy
                cur_idx += 1
            return cur_idx

        def find_min_sym_idx() -> int:
            """Find the minimum index of the symmetry type

            Returns:
                int: Minimum index of the symmetry type
            """

            def bisect_l(data: List[DataMO], key: str) -> int:
                min_idx, max_idx = 0, len(data)
                while min_idx < max_idx:
                    mid_idx = (min_idx + max_idx) // 2
                    if data[mid_idx].sym_type < key:
                        min_idx = mid_idx + 1
                    else:
                        max_idx = mid_idx
                return min_idx

            key = symmetry_type
            min_idx = bisect_l(data_all_mo.electronic, key)
            return min_idx

        def create_within_moltra_dict(start_idx: int, end_idx: int) -> OrderedDict:
            """Create the within_moltra dictionary

            Returns:
                OrderedDict: within_moltra dictionary
            """
            within_moltra = OrderedDict(sorted({mo.eigenvalue_no: False for mo in data_all_mo.electronic}.items()))
            for idx in range(start_idx, end_idx):
                within_moltra[data_all_mo.electronic[idx].eigenvalue_no] = True
            return within_moltra

        def create_energy_str(within_moltra: ODict[int, bool]):
            energy_str = ""
            cur_mo_num = 0
            start_mo_num = cur_mo_num
            left = True
            ser_end = False
            for eigenvalue_no, is_used in within_moltra.items():
                cur_mo_num = eigenvalue_no
                if is_used:
                    if left:
                        if energy_str == "":
                            energy_str += f"{cur_mo_num}"
                        else:
                            energy_str += f",{cur_mo_num}"
                        start_mo_num = cur_mo_num
                        left = False
                        ser_end = True
                elif ser_end:
                    if cur_mo_num > start_mo_num + 1:
                        energy_str += f"..{cur_mo_num - 1}"
                    left = True
                    ser_end = False
            if ser_end:
                energy_str += f"..{cur_mo_num}"
            return energy_str

        start_sym_idx = find_min_sym_idx()
        energy_str = energy_str.upper().replace("ENERGY", "")
        min_energy, max_energy, step = map(float, energy_str.split())
        if min_energy > max_energy:
            msg = f"The minimum energy is larger than the maximum energy: {min_energy} > {max_energy}"
            raise ValueError(msg)
        min_energy_idx = get_min_energy_idx(min_energy, step)
        max_energy_idx = get_max_energy_idx(max_energy, step)

        within_moltra = create_within_moltra_dict(min_energy_idx + start_sym_idx, max_energy_idx + start_sym_idx)
        energy_str = create_energy_str(within_moltra)
        return energy_str

    def __parse_range_str(self, range_str: str, symmetry_type: str) -> str:
        """Parse the range string
        (e.g.) "10..180, 200..300, 400..oo" => "10..180,200..300,400..500"

        Args:
            range_str (str): Range string

        Returns:
            str: Range string
        """
        # [-]?(?:[0-9]+|oo) => -oo or oo or integer
        # \.{2} => ..
        regex = r"[-]?(?:[0-9]+|oo)\.{2}[-]?(?:[0-9]+|oo)"
        items: List[str] = re.findall(regex, range_str)
        ret_list: List[str] = []
        for item in items:
            item_li = item.replace("..", " ").split()
            range_li = [1 if item == "-oo" else len(self.eigenvalues.energies[symmetry_type]) if item == "oo" else int(item) for item in item_li]
            if range_li[0] > range_li[1]:
                msg = f"The minimum index is larger than the maximum index: {range_li[0]} > {range_li[1]}\n\
your input: {range_str}, invalid input part: {item}"
                raise ValueError(msg)
            ret_list.append(f"{range_li[0]}..{range_li[1]}")
        return ",".join(ret_list)
