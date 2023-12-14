import re
from bisect import bisect_left, bisect_right
from io import TextIOWrapper
from typing import List

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
        """
        dirac_output.seek(0)
        self.__read_electron_number(dirac_output)
        dirac_output.seek(0)
        self.__read_moltra(dirac_output)
        self.__read_eigenvalues(dirac_output)
        self.__duplicate_moltra_str()
        self.__calculate_moltra_idx_range()

    def __read_electron_number(self, dirac_output: TextIOWrapper) -> None:
        self.electrons = get_electron_num_from_input(dirac_output)
        if self.electrons == 0:
            self.electrons = get_electron_num_from_scf_field(dirac_output)

    def __read_eigenvalues(self, dirac_output: TextIOWrapper) -> None:
        self.eigenvalues.get_eigenvalues(dirac_output)

    def __read_moltra(self, dirac_output: TextIOWrapper) -> None:
        self.moltra_info.read_moltra_section(dirac_output)

    def __duplicate_moltra_str(self) -> None:
        """Duplicate the moltra range string if it is not enough"""
        if self.moltra_info.is_default:
            # Set the default range string
            for _ in range(len(self.eigenvalues.shell_num)):
                self.moltra_info.range_str.append("ENERGY -20.0 10.0 1.0")
        if len(self.moltra_info.range_str) != len(self.eigenvalues.shell_num):
            # one-line input is allowed, duplicate the first line
            # https://gitlab.com/dirac/dirac/-/blob/ea717cdb294035d8af3ebe2b1e00cf94f1c1a6b7/src/moltra/trainp.F#L592-600
            for _ in range(len(self.eigenvalues.shell_num) - len(self.moltra_info.range_str)):
                self.moltra_info.range_str.append(self.moltra_info.range_str[0])

    def __calculate_moltra_idx_range(self) -> None:
        keys = list(self.eigenvalues.shell_num.keys())
        for i, item in enumerate(self.moltra_info.range_str):
            symmetry_type = keys[i]
            if "ALL" == item.upper():
                self.moltra_info.range_str[i] = f"1..{len(self.eigenvalues.energies[symmetry_type])}"
            elif "ENERGY" in item.upper():
                self.moltra_info.range_str[i] = self.__parse_energy_str(item, symmetry_type)
            else:
                self.moltra_info.range_str[i] = self.__parse_range_str(item, symmetry_type)

    def __parse_energy_str(self, energy_str: str, symmetry_type: str) -> str:
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

        energy_str = energy_str.replace("ENERGY", "")
        min_energy, max_energy, step = map(float, energy_str.split())
        if min_energy > max_energy:
            msg = f"The minimum energy is larger than the maximum energy: {min_energy} > {max_energy}"
            raise ValueError(msg)

        min_energy_idx = get_min_energy_idx(min_energy, step)
        max_energy_idx = get_max_energy_idx(max_energy, step)
        return f"{min_energy_idx+1}..{max_energy_idx}"

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
        tmp_list: List[str] = []
        for item in items:
            item_li = item.replace("..", " ").split()
            range_li = [1 if item == "-oo" else len(self.eigenvalues.energies[symmetry_type]) if item == "oo" else int(item) for item in item_li]
            if range_li[0] > range_li[1]:
                msg = f"The minimum index is larger than the maximum index: {range_li[0]} > {range_li[1]}\n\
your input: {range_str}, invalid input part: {item}"
                raise ValueError(msg)
            tmp_list.append(f"{range_li[0]}..{range_li[1]}")
        return ",".join(tmp_list)
