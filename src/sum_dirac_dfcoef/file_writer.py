from pathlib import Path
from typing import List

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.data import DataMO
from sum_dirac_dfcoef.header_info import HeaderInfo
from sum_dirac_dfcoef.utils import debug_print


class OutputFileWriter:
    def __init__(self) -> None:
        super().__init__()
        self.output_path = self.get_output_path()

    def write_headerinfo(self, header_info: HeaderInfo) -> None:
        with open(self.output_path, "a", encoding="utf-8") as f:
            line = ""
            line += f"electron_num {header_info.electrons} "
            for symmetry_type, d in header_info.eigenvalues.shell_num.items():
                line += f"{symmetry_type} "
                for eigenvalue_type, num in d.items():
                    line += f"{eigenvalue_type} {num} "
            line += "\n"
            f.write(line)

    def write_mo_data(self, mo_data: List[DataMO], add_blank_line: bool) -> None:
        with open(self.output_path, "a", encoding="utf-8") as f:
            for mo in mo_data:
                digit_int = len(str(int(mo.mo_energy)))  # number of digits of integer part
                # File write but if args.compress is True \n is not added
                mo_info_energy = f"{mo.mo_info} {mo.mo_energy:{digit_int}.{args.decimal}f}" + ("\n" if not args.compress else "")
                f.write(mo_info_energy)

                for c in mo.coef_list:
                    for idx in range(c.multiplication):
                        percentage = c.coefficient / mo.norm_const_sum * 100
                        atomic_symmetry_label = f"{c.function_label}({c.start_idx + idx})" if c.need_identifier else c.function_label
                        output_str: str
                        if args.compress:
                            output_str = f" {atomic_symmetry_label} {percentage:.{args.decimal}f}"
                        else:
                            output_str = f"{atomic_symmetry_label:<12} {percentage:{args.decimal+4}.{args.decimal}f} %\n"
                        f.write(output_str)
                f.write("\n")  # add empty line
                debug_print(f"sum of coefficient {mo.norm_const_sum:.{args.decimal}f}")
            if add_blank_line:
                f.write("\n")

    def create_blank_file(self) -> None:
        # Open the file in write mode
        # Even if the file already exists, it will be overwritten with a blank file
        file = open(self.output_path, "w", encoding="utf-8")
        file.close()

    def get_output_path(self) -> Path:
        if args.output is None:
            output_name = "sum_dirac_dfcoef.out"
            output_path = Path.absolute(Path.cwd() / output_name)
        else:
            output_name = args.output
            if Path(output_name).is_absolute():
                output_path = Path(output_name)
            else:
                output_path = Path.absolute(Path.cwd() / output_name)
        return output_path


output_file_writer = OutputFileWriter()
