#!/usr/bin/env python3

import copy
import sys
from pathlib import Path
from typing import List

from sum_dirac_dfcoef.args import args
from sum_dirac_dfcoef.atoms import AtomInfo
from sum_dirac_dfcoef.coefficient import get_coefficient
from sum_dirac_dfcoef.data import DataAllMO, DataMO
from sum_dirac_dfcoef.file_writer import output_file_writer
from sum_dirac_dfcoef.functions_info import get_functions_info
from sum_dirac_dfcoef.header_info import HeaderInfo
from sum_dirac_dfcoef.utils import debug_print, space_separated_parsing


def is_this_row_for_coefficients(words: List[str]) -> bool:
    # min: 4 coefficients and other words => 5 words
    if 5 <= len(words) <= 9 and words[0].isdigit():
        return True
    else:
        return False


def need_to_skip_this_line(words: List[str]) -> bool:
    if len(words) <= 1:
        return True
    else:
        return False


def need_to_create_results_for_current_mo(words: List[str], is_reading_coefficients: bool) -> bool:
    if is_reading_coefficients and len(words) <= 1:
        return True
    else:
        return False


def need_to_get_mo_sym_type(words: List[str], start_mo_coefficients: bool) -> bool:
    if not start_mo_coefficients and len(words) == 3 and words[0] == "Fermion" and words[1] == "ircop":
        return True
    return False


def need_to_start_mo_section(words: List[str], start_mo_coefficients: bool) -> bool:
    if not start_mo_coefficients and words[1] == "Electronic" and words[2] == "eigenvalue" and "no." in words[3]:
        return True
    elif not start_mo_coefficients and words[1] == "Positronic" and words[2] == "eigenvalue" and "no." in words[3]:
        return True
    return False


def get_dirac_filepath() -> Path:
    if not args.file:
        sys.exit("ERROR: DIRAC output file is not given. Please use -f option.")
    elif not Path(args.file).exists():
        sys.exit(f"ERROR: DIRAC output file is not found. file={args.file}")
    path = Path(args.file)
    return path


def check_start_vector_print(words: List[str]) -> bool:
    # ****************************** Vector print ******************************
    if len(words) < 4:
        return False
    elif words[1] == "Vector" and words[2] == "print":
        return True
    return False


def check_end_vector_print(
    words: List[str],
    start_vector_print: bool,
    start_mo_section: bool,
    start_mo_coefficients: bool,
    is_reading_coefficients: bool,
) -> bool:
    # https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/issues/7#issuecomment-1377969626
    if len(words) >= 2 and start_vector_print and start_mo_section and not start_mo_coefficients and not is_reading_coefficients:
        return True
    return False


def should_write_positronic_results_to_file() -> bool:
    if args.all_write or args.positronic_write:
        return True
    else:
        return False


def should_write_electronic_results_to_file() -> bool:
    if args.all_write or not args.positronic_write:
        return True
    else:
        return False


def main() -> None:
    is_electronic: bool = False
    is_reading_coefficients: bool = False
    start_mo_coefficients: bool = False
    start_mo_section: bool = False
    start_vector_print: bool = False
    mo_sym_type: str = ""

    dirac_filepath = get_dirac_filepath()
    dirac_output = open(dirac_filepath, encoding="utf-8")
    dirac_output.seek(0)  # rewind to the beginning of the file
    if (not args.no_scf):
        header_info = HeaderInfo()
        header_info.read_header_info(dirac_output)
    dirac_output.seek(0)
    functions_info = get_functions_info(dirac_output)
    output_file_writer.create_blank_file()
    if (args.no_scf or args.positronic_write):
        # If args.no_scf is True, we cannot get header_info from the output file of DIRAC.
        # also, if args.positronic_write is True, we do not need to write header_info to the output file.
        # because dcaspt2_input_generator program does not support positronic results.
        output_file_writer.write_no_header_info()
    else:
        output_file_writer.write_headerinfo(header_info)

    data_mo = DataMO()
    data_all_mo = DataAllMO()
    used_atom_info: dict[str, AtomInfo] = {}
    current_atom_info = AtomInfo()
    dirac_output.seek(0)  # rewind to the beginning of the file
    for line_str in dirac_output:
        words: List[str] = space_separated_parsing(line_str)

        if not start_vector_print:
            if check_start_vector_print(words):
                start_vector_print = True
            continue

        if need_to_get_mo_sym_type(words, start_mo_coefficients):
            mo_sym_type = words[2]

        elif need_to_skip_this_line(words):
            if need_to_create_results_for_current_mo(words, is_reading_coefficients):
                start_mo_coefficients = False
                data_mo.fileter_coefficients_by_threshold()
                if is_electronic:
                    data_all_mo.electronic.append(copy.deepcopy(data_mo))
                else:  # Positronic
                    data_all_mo.positronic.append(copy.deepcopy(data_mo))
                debug_print(f"End of reading {data_mo.electron_num}th MO")
                is_reading_coefficients = False

        elif need_to_start_mo_section(words, start_mo_coefficients):
            """
            (e.g.)
            words = ["*", "Electronic", "eigenvalue", "no.", "22:", "-2.8417809384721"]
            words = ["*", "Electronic", "eigenvalue", "no.122:", "-2.8417809384721"]
            """
            start_mo_section = True
            start_mo_coefficients = True
            if words[1] == "Positronic":
                is_electronic = False
            elif words[1] == "Electronic":
                is_electronic = True
            else:
                msg = f"UnKnow MO type, MO_Type={words[1]}"
                raise Exception(msg)
            try:
                electron_num = int(words[-2][:-1].replace("no.", ""))
            except ValueError:
                # If *** is printed, we have no information about what number this MO is.
                # Therefore, we assume that electron_num is the next number after prev_electron_num.
                prev_electron_num = data_mo.electron_num  # prev_electron is the number of electrons of the previous MO
                electron_num = prev_electron_num + 1
            mo_energy = float(words[-1])
            mo_info = (
                f"{mo_sym_type} {electron_num}"
                if args.compress
                else (f"Electronic no. {electron_num} {mo_sym_type}" if is_electronic else f"Positronic no. {electron_num} {mo_sym_type}")
            )
            # Here is the start point of reading coefficients of the current MO
            data_mo.reset()  # reset data_mo because we need to delete data_mo of the previous MO
            data_mo.electron_num = electron_num
            data_mo.mo_energy = mo_energy
            data_mo.mo_info = mo_info
            used_atom_info.clear()  # reset used_atom_info because we need to delete used_atom_info of the previous MO

        elif check_end_vector_print(
            words,
            start_vector_print,
            start_mo_section,
            start_mo_coefficients,
            is_reading_coefficients,
        ):
            # End of reading coefficients
            break

        # Read coefficients or the end of coefficients section
        elif start_mo_coefficients:
            if not is_this_row_for_coefficients(words):
                continue
            is_reading_coefficients = True
            component_func = "large" if line_str[10] == "L" else ("small" if line_str[10] == "S" else "")  # CLS
            symmetry_label = line_str[12:15].strip()  # REP (e.g. "Ag "), symmetry_label="Ag"
            atom_label = line_str[15:18].strip()  # NAMN (e.g. "Cm "), atom_labe="Cm"
            gto_type = line_str[18:22].strip()  # GTOTYP (e.g. "s   "), gto_type="s"
            label = symmetry_label + atom_label

            if current_atom_info.count_remaining_functions() == 0 or label != current_atom_info.label:
                # First, we need to read information about the current atom.
                if label not in used_atom_info:
                    # It is the first time to read information about the current atom.
                    cur_atom_start_idx = 1
                else:
                    # It is not the first time to read information about the current atom.
                    # So we need to read information about the previous atom from used_atom_info.
                    # current start_idx = previous start_idx + previous mul
                    cur_atom_start_idx = used_atom_info[label].start_idx + used_atom_info[label].mul
                # Validate start_idx
                if cur_atom_start_idx not in functions_info[component_func][symmetry_label][atom_label]:
                    msg = f"start_idx={cur_atom_start_idx} is not found in functions_info[{component_func}][{symmetry_label}][{atom_label}]"
                    raise Exception(msg)
                # We can get information about the current atom from functions_info with start_idx.
                current_atom_info = copy.deepcopy(functions_info[component_func][symmetry_label][atom_label][cur_atom_start_idx])
                # Update used_atom_info with current_atom_info
                used_atom_info[label] = copy.deepcopy(current_atom_info)

            current_atom_info.decrement_function(gto_type)
            data_mo.add_coefficient(get_coefficient(line_str, functions_info, current_atom_info.start_idx))

    # End of reading file
    # Write results to the file
    if not args.no_sort:
        data_all_mo.electronic.sort(key=lambda x: x.mo_energy)
        data_all_mo.positronic.sort(key=lambda x: x.mo_energy)
    if should_write_positronic_results_to_file():  # Positronic
        # Write positronic results to the file
        output_file_writer.write_mo_data(data_all_mo.positronic, add_blank_line=True)
    if should_write_electronic_results_to_file():  # Electronic
        # Write electronic results to the file
        output_file_writer.write_mo_data(data_all_mo.electronic, add_blank_line=False)
