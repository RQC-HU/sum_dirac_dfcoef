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
    if not args.no_scf:
        header_info = HeaderInfo()
        header_info.read_header_info(dirac_output)
    dirac_output.seek(0)
    functions_info = get_functions_info(dirac_output)
    output_file_writer.create_blank_file()
    if args.no_scf or args.positronic_write:
        # If args.no_scf is True, we cannot get header_info from the output file of DIRAC.
        # also, if args.positronic_write is True, we do not need to write header_info to the output file
        # because dcaspt2_input_generator program does not support positronic results.
        output_file_writer.write_no_header_info()
    else:
        output_file_writer.write_headerinfo(header_info)

    # Read coefficients from the output file of DIRAC and store them in data_all_mo.
    privec_processor = PrivecProcessor(dirac_output, functions_info)
    privec_processor.read_privec_data()
    data_all_mo = privec_processor.data_all_mo

    if should_write_electronic_results_to_file():
        add_blank_line = True if args.all_write else False
        output_file_writer.write_mo_data(data_all_mo.electronic, add_blank_line=add_blank_line)
    if should_write_positronic_results_to_file():
        output_file_writer.write_mo_data(data_all_mo.positronic, add_blank_line=False)
