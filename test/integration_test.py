import os
import re
import subprocess
from pathlib import Path
from typing import List

import pytest


class Env:
    def __init__(self, input_filename: str, options: str, ref_filename: str = "", result_filename: str = "") -> None:
        self.test_path = Path(__file__).resolve().parent
        self.ref_filepath = Path.joinpath(self.test_path, "references", ref_filename)
        self.result_filepath = Path.joinpath(self.test_path, "results", result_filename)
        self.input_filepath = Path.joinpath(self.test_path, "data", input_filename)
        if ref_filename == "":
            self.command: str = f"sum_dirac_dfcoef -i {self.input_filepath} {options}"
        else:
            self.command: str = f"sum_dirac_dfcoef -i {self.input_filepath} -o {self.result_filepath} {options}"


def get_output_list(filepath: Path) -> List[List[str]]:
    return [re.split(" +", line.rstrip("\r\n")) for line in list(open(filepath).read().splitlines())]


@pytest.mark.parametrize(
    "ref_filename, result_filename, input_filename, options",
    # fmt: off
    [
        ("ref.Ar.compress.out"                      , "result.Ar.compress.out"                      , "Ar_Ar.out"                    , "-d 15 -c"),
        ("ref.Ar.no_sort.compress.out"              , "result.Ar.no_sort.compress.out"              , "Ar_Ar.out"                    , "-d 15 --no-sort -c"),
        ("ref.N2.compress.out"                      , "result.N2.compress.out"                      , "N2_N2.out"                    , "-d 15 -c"),
        ("ref.N2.compress.positronic.out"           , "result.N2.compress.positronic.out"           , "N2_N2.out"                    , "-d 15 -pc"),
        ("ref.N2.compress.all.out"                  , "result.N2.compress.all.out"                  , "N2_N2.out"                    , "-d 15 -ac"),
        ("ref.N2.no_sort.compress.out"              , "result.N2.no_sort.compress.out"              , "N2_N2.out"                    , "-d 15 --no-sort -c"),
        ("ref.N2.no_sort.compress.positronic.out"   , "result.N2.no_sort.compress.positronic.out"   , "N2_N2.out"                    , "-d 15 --no-sort -pc"),
        ("ref.N2.no_sort.compress.all.out"          , "result.N2.no_sort.compress.all.out"          , "N2_N2.out"                    , "-d 15 --no-sort -ac"),
        ("ref.uo2.no_sort.compress.out"             , "result.uo2.no_sort.compress.out"             , "x2c_uo2_238.out"              , "-d 15 --no-sort -c"),
        ("ref.ucl4.compress.out"                    , "result.ucl4.compress.out"                    , "x2c_ucl4.out"                 , "-d 15 -g"),
        ("ref.ucl4.no_sort.compress.out"            , "result.ucl4.no_sort.compress.out"            , "x2c_ucl4.out"                 , "-d 15 --no-sort -c"),
        ("ref.Cm3+_phen.compress.out"               , "result.Cm3+_phen.compress.out"               , "x2c_Cm3+_phen.out"            , "-d 15 -g"),
        ("ref.Cm3+_phen_reorder.compress.out"       , "result.Cm3+_phen_reorder.compress.out"       , "x2c_Cm3+_phen_reorder.out"    , "-d 15 -g"),
        ("ref.H2.no_scf.compress.out"               , "result.H2.no_scf.compress.out"               , "H2.noscf_H2.out"              , "-d 15 -c --no-scf"),
        ("ref.H2O.invalid.eigpri.compress.out"      , "result.H2O.invalid.eigpri.compress.out"      , "H2O.invalid.eigpri.out"       , "-d 15 -c --no-scf"),
        ("ref.methane.whitespace.compress.out"      , "result.methane.whitespace.compress.out"      , "methane.whitespace_mol.out"   , "-d 15 -c"),
        # multiprocess (should be the same as the single process case)
        ("ref.ucl4.compress.out"                    , "result.ucl4.compress.multi-process.out"      , "x2c_ucl4.out"                 , "-j2 -d 15 -g"),
        # DIRAC 19 UO2 x2c (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/93#issue-2164290127)
        ("ref.uo2.compress.out"                     , "result.uo2.compress.out"                     , "x2c_uo2_238.out"              , "-d 15 -g"),
        # DIRAC 19 N2 4component (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/93#issuecomment-1976947902)
        ("ref.N2.dirac19.compress.out"              , "result.N2.dirac19.compress.out"              , "dirac19_N2_N2.out"            , "-d 15 -g"),
        # F2 (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/92)
        ("ref.F2.compress.out"                      , "result.F2.compress.out"                      , "F2_F2.out"                    , "-d 15 -g"),
        # CO (No inversion symmetry linear molecule, https://github.com/RQC-HU/sum_dirac_dfcoef/issues/95)
        ("ref.CO.compress.out"                      , "result.CO.compress.out"                      , "CO_CO.out"                    , "-d 15 -g"),
        # Explicitly set **MOLTRA > .SCHEME option (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/108)
        ("ref.N2.scheme4.dirac22.compress.out"      , "result.N2.scheme4.dirac22.compress.out"      , "dirac22_scheme4_N2.out"       , "-d 15 -g"),
        ("ref.N2.scheme6.dirac23.compress.out"      , "result.N2.scheme6.dirac23.compress.out"      , "dirac23_scheme6_N2.out"       , "-d 15 -g"),
        # Ignore magnetic quantum number label (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/115)
        ("ref.ucl4.ignore_ml.compress.out"          , "result.ucl4.ignore_ml.compress.out"          , "x2c_ucl4.out"                 , "-d 15 -g --ignore-ml"),
        # Ignore symmetry label (https://github.com/RQC-HU/sum_dirac_dfcoef/issues/117)
        ("ref.uo2.ignore_sym.compress.out"          , "result.uo2.ignore_sym.compress.out"          , "x2c_uo2_238.out"              , "-d 15 -c --ignore-sym"),
        # Both ignore labels
        ("ref.ucl4.ignore_ml_and_sym.compress.out"  , "result.ucl4.ignore_ml_and_sym.compress.out"  , "x2c_ucl4.out"                 , "-d 15 -g --ignore-ml --ignore-sym"),
    ]
    # fmt: on
)
def test_sum_dirac_dfcoeff_compress(ref_filename: str, result_filename: str, input_filename: str, options: str):
    env = Env(input_filename, options, ref_filename, result_filename)
    os.chdir(env.test_path)
    print(f"{env.test_path} test start...\ncommand: {env.command}")
    subprocess.run(env.command.split(), encoding="utf-8", check=True)

    ref_list: List[List[str]] = get_output_list(env.ref_filepath)
    result_list: List[List[str]] = get_output_list(env.result_filepath)
    # File should have the same number of lines
    assert len(ref_list) == len(result_list), f"Number of lines in {ref_filename}(={len(ref_list)}) and {result_filename}(={len(result_list)}) are different."
    threshold: float = 1e-10
    header = True
    for line_idx, (ref, out) in enumerate(zip(ref_list, result_list)):
        if header:
            assert ref == out
            if len(ref) == 1 and len(out) == 1:  # End header
                header = False
        else:
            if len(ref) < 2 or len(out) < 2:
                continue
            # ref[0]: irrep, ref[1]: energy order index in the irrep, ref[2]: energy, ref[3:]: Symmetry value and coefficient
            # (e.g.) E1u 19 -8.8824415703374 B3uUpx 49.999172476298732 B2uUpy 49.999172476298732
            assert ref[0] == out[0], f"irrep in line {line_idx} of {ref_filename} and {result_filename} are different."
            assert ref[1] == out[1], f"Energy order index in line {line_idx} of {ref_filename} and {result_filename} are different."
            assert float(ref[2]) == pytest.approx(float(out[2]), abs=threshold), f"Energy in line {line_idx} of {ref_filename} and {result_filename} are different."
            for idx, (ref_val, out_val) in enumerate(zip(ref[3:], out[3:])):
                if idx % 2 == 0:
                    assert ref_val == out_val, f"Symmetry value in line {line_idx} of {ref_filename} and {result_filename} are different."
                else:
                    assert float(ref_val) == pytest.approx(
                        float(out_val), abs=threshold
                    ), f"Contribution of the AO in the MO in line {line_idx} of {ref_filename} and {result_filename} are different."


@pytest.mark.parametrize(
    "ref_filename, result_filename, input_filename, options",
    # fmt: off
    [
        ("ref.Ar.out"                       , "result.Ar.out"                       , "Ar_Ar.out"                       , "-d 15"),
        ("ref.Ar.no_sort.out"               , "result.Ar.no_sort.out"               , "Ar_Ar.out"                       , "-d 15 --no-sort"),
        ("ref.N2.out"                       , "result.N2.out"                       , "N2_N2.out"                       , "-d 15"),
        ("ref.N2.positronic.out"            , "result.N2.positronic.out"            , "N2_N2.out"                       , "-d 15 -p"),
        ("ref.N2.all.out"                   , "result.N2.all.out"                   , "N2_N2.out"                       , "-d 15 -a"),
        ("ref.N2.no_sort.out"               , "result.N2.no_sort.out"               , "N2_N2.out"                       , "-d 15 --no-sort"),
        ("ref.N2.no_sort.positronic.out"    , "result.N2.no_sort.positronic.out"    , "N2_N2.out"                       , "-d 15 --no-sort -p"),
        ("ref.N2.no_sort.all.out"           , "result.N2.no_sort.all.out"           , "N2_N2.out"                       , "-d 15 --no-sort -a"),
        ("ref.uo2.special.out"              , "result.uo2.special.out"              , "special_exit_condition_UO2.out"  , "-d 15"),
        ("ref.uo2.out"                      , "result.uo2.out"                      , "x2c_uo2_238.out"                 , "-d 15"),
        ("ref.uo2.no_sort.out"              , "result.uo2.no_sort.out"              , "x2c_uo2_238.out"                 , "-d 15 --no-sort"),
        ("ref.ucl4.out"                     , "result.ucl4.out"                     , "x2c_ucl4.out"                    , "-d 15"),
        ("ref.ucl4.no_sort.out"             , "result.ucl4.no_sort.out"             , "x2c_ucl4.out"                    , "-d 15 --no-sort"),
        ("ref.Cm3+_phen.out"                , "result.Cm3+_phen.out"                , "x2c_Cm3+_phen.out"               , "-d 15"),
        # multiprocess (should be the same as the single process case)
        ("ref.uo2.out"                      , "result.uo2.multi-process.out"        , "x2c_uo2_238.out"                 , "-j2 -d 15"),
    ]
    # fmt: on
)
def test_sum_dirac_dfcoeff(ref_filename: str, result_filename: str, input_filename: str, options: str):
    env = Env(input_filename, options, ref_filename, result_filename)
    os.chdir(env.test_path)
    print(f"{env.test_path} test start...\ncommand: {env.command}")
    subprocess.run(env.command.split(), encoding="utf-8", check=True)

    ref_list: List[List[str]] = get_output_list(env.ref_filepath)
    result_list: List[List[str]] = get_output_list(env.result_filepath)
    # File should have the same number of lines
    assert len(ref_list) == len(result_list), f"Number of lines in {ref_filename}(={len(ref_list)}) and {result_filename}(={len(result_list)}) are different."
    threshold: float = 1e-10
    header = True
    for line_idx, (ref, out) in enumerate(zip(ref_list, result_list)):
        if header:
            assert ref == out
            if len(ref) == 1 and len(out) == 1:  # End header
                header = False
        else:
            if len(ref) < 2 or len(out) < 2:
                continue
            if "%" in ref[-1]:
                ref_value = float(ref[-2])
                out_value = float(out[-2])
                ref_list_str = " ".join(ref[:-2])
                out_list_str = " ".join(out[:-2])
            else:
                ref_value = float(ref[-1])
                out_value = float(out[-1])
                ref_list_str = " ".join(ref[:-1])
                out_list_str = " ".join(out[:-1])
            assert ref_list_str == out_list_str, f"line {line_idx}: {ref_list_str} != {out_list_str}\nref: {ref_list[line_idx]}\nout:{result_list[line_idx]}"
            assert ref_value == pytest.approx(out_value, abs=threshold), f"line {line_idx}: {ref_value} != {out_value}\nref: {ref_list[line_idx]}\nout:{result_list[line_idx]}"


@pytest.mark.parametrize(
    "input_filename, options, expected_error_message",
    # fmt: off
    [
        ("H2O.invalid.eigpri.out"  , "-g -d 15"    , "Your .EIGPRI option in your DIRAC input file is invalid!"),
        ("H2.noscf_H2.out"         , "-g -d 15"    , "Cannot find SCF calculation settings"),
        ("Ar_Ar.out"               , "-g --no-scf" , "-g/--for-generator and --no-scf options cannot be set at the same time"),
        ("Ar_Ar.out"               , "-g -p"       , "-g/--for-generator and -p/--positronic-write options cannot be set at the same time"),
    ],
    # fmt: on
)
def test_invalid_option_raise_error(input_filename: str, options: str, expected_error_message: str):
    env = Env(input_filename, options)
    os.chdir(env.test_path)
    print(f"{env.test_path} test start...\ncommand: {env.command}")
    # Run the command and capture the output
    p = subprocess.run(env.command.split(), encoding="utf-8", stderr=subprocess.PIPE, check=False)
    assert p.returncode != 0, "Failed: Command should fail but return code is 0."
    assert expected_error_message in p.stderr, f"Failed: Expected error message not found in output:\nexptected: {expected_error_message}\nstderr: {p.stderr}"


@pytest.mark.parametrize(
    "ref_filename, result_filename, input_filename, options, expected_warning_message",
    # fmt: off
    [
        ("ref.Ar.no_vector_print_data.out", "result.Ar.no_vector_print_data.out", "no_vector_print_data_Ar_Ar.out", "-d 15 -c", "WARNING: The next title is detected before the end of reading coefficients."),  # noqa: E501
    ],
    # fmt: on
)
def test_no_vector_print_data(ref_filename: str, result_filename: str, input_filename: str, options: str, expected_warning_message: str):
    env = Env(input_filename, options, ref_filename, result_filename)
    os.chdir(env.test_path)
    print(f"{env.test_path} test start...\ncommand: {env.command}")
    process = subprocess.run(env.command.split(), encoding="utf-8", check=True, stdout=subprocess.PIPE)

    # Check the warning message
    assert expected_warning_message in process.stdout, f"Failed: Expected warning message not found in output:\nexptected: {expected_warning_message}\nstdout: {process.stdout}"

    # Check the output file
    ref_list: List[List[str]] = get_output_list(env.ref_filepath)
    result_list: List[List[str]] = get_output_list(env.result_filepath)
    # File should have the same number of lines
    assert len(ref_list) == len(result_list), f"Number of lines in {ref_filename}(={len(ref_list)}) and {result_filename}(={len(result_list)}) are different."
    threshold: float = 1e-10
    header = True
    for line_idx, (ref, out) in enumerate(zip(ref_list, result_list)):
        if header:
            assert ref == out
            if len(ref) == 1 and len(out) == 1:  # End header
                header = False
        else:
            if len(ref) < 2 or len(out) < 2:
                continue
            # ref[0]: irrep, ref[1]: energy order index in the irrep, ref[2]: energy, ref[3:]: Symmetry value and coefficient
            # (e.g.) E1u 19 -8.8824415703374 B3uUpx 49.999172476298732 B2uUpy 49.999172476298732
            assert ref[0] == out[0], f"irrep in line {line_idx} of {ref_filename} and {result_filename} are different."
            assert ref[1] == out[1], f"Energy order index in line {line_idx} of {ref_filename} and {result_filename} are different."
            assert float(ref[2]) == pytest.approx(float(out[2]), abs=threshold), f"Energy in line {line_idx} of {ref_filename} and {result_filename} are different."
            for idx, (ref_val, out_val) in enumerate(zip(ref[3:], out[3:])):
                if idx % 2 == 0:
                    assert ref_val == out_val, f"Symmetry value in line {line_idx} of {ref_filename} and {result_filename} are different."
                else:
                    assert float(ref_val) == pytest.approx(
                        float(out_val), abs=threshold
                    ), f"Contribution of the AO in the MO in line {line_idx} of {ref_filename} and {result_filename} are different."


def test_version_option():
    command = "sum_dirac_dfcoef -v"
    p = subprocess.run(command.split(), encoding="utf-8", check=True, stdout=subprocess.PIPE)
    # p.stdout = x.y.z\n
    out_version_str = p.stdout.rstrip("\n")
    about_file = Path(__file__).resolve().parent.parent / "src/sum_dirac_dfcoef/__about__.py"
    # readline() => __version__ = "x.y.z"
    # split()[-1] => "x.y.z"
    # replace() => x.y.z
    ref_version_str = open(about_file).readline().split()[-1].replace('"', "")
    assert out_version_str == ref_version_str


def test_help_option():
    command = "sum_dirac_dfcoef -h"
    p = subprocess.run(command.split(), check=True)
    assert p.returncode == 0
