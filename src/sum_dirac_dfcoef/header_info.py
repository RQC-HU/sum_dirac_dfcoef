from io import TextIOWrapper

from sum_dirac_dfcoef.eigenvalues import Eigenvalues
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

    def read_header_info(self, dirac_output: TextIOWrapper) -> None:
        """Read the header information from the output file of DIRAC

        Args:
            dirac_output (TextIOWrapper): Output file of DIRAC

        Returns:
        """
        self.__read_eigenvalues(dirac_output)
        self.__read_moltra(dirac_output)
        self.__validate_len_moltra()

    def __read_eigenvalues(self, dirac_output: TextIOWrapper) -> None:
        self.eigenvalues.get_eigenvalues(dirac_output)

    def __read_moltra(self, dirac_output: TextIOWrapper) -> None:
        self.moltra_info.read_moltra_section(dirac_output)

    def __validate_len_moltra(self) -> None:
        if self.moltra_info.is_default:
            return
        if len(self.moltra_info.ranges) != len(self.eigenvalues.shell_num):
            msg = f"The number of lines in the MOLTRA .ACTIVE section is not equal to the number of the symmetry types\
in the Eigenvalues section.\n\
symmetry types in the Eigenvalues section: {self.eigenvalues.shell_num.keys()}\n\
lines in the MOLTRA .ACTIVE section: {self.moltra_info.ranges}"
            raise ValueError(msg)
