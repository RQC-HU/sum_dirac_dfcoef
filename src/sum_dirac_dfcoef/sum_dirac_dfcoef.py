#!/usr/bin/env python3
from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.file_writer import output_file_writer
from sum_dirac_dfcoef.functions_info import get_functions_info
from sum_dirac_dfcoef.header_info import HeaderInfo
from sum_dirac_dfcoef.privec_reader import PrivecProcessor
from sum_dirac_dfcoef.utils import get_dirac_filepath, should_write_electronic_results_to_file, should_write_positronic_results_to_file


def main() -> None:
    dirac_filepath = get_dirac_filepath()
    dirac_output = open(dirac_filepath, encoding="utf-8")
    dirac_output.seek(0)  # rewind to the beginning of the file
    header_info = HeaderInfo()
    if args.for_generator:
        header_info.read_header_info(dirac_output)
    dirac_output.seek(0)
    functions_info = get_functions_info(dirac_output)
    output_file_writer.create_blank_file()
    if not args.for_generator:
        # If the output file is not for dcaspt2_input_generator, don't write header information.
        output_file_writer.write_no_header_info()
    else:
        output_file_writer.write_headerinfo(header_info)

    # Read coefficients from the output file of DIRAC and store them in data_all_mo.
    privec_processor = PrivecProcessor(dirac_output, functions_info, header_info.eigenvalues)
    privec_processor.read_privec_data()
    privec_processor.data_all_mo.sort_mo_energy()
    data_all_mo = privec_processor.data_all_mo

    if should_write_electronic_results_to_file():
        add_blank_line = True if args.all_write else False
        output_file_writer.write_mo_data(data_all_mo.electronic, add_blank_line=add_blank_line)
    if should_write_positronic_results_to_file():
        output_file_writer.write_mo_data(data_all_mo.positronic, add_blank_line=False)
