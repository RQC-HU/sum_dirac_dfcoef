import os
import re
import subprocess
from typing import List

import pytest


@pytest.mark.parametrize(
    "ref_filename, result_filename, input_filename, options",
    # fmt: off
    [
        ("ref.Ar.compress.out"                      , "result.Ar.compress.out"                      , "Ar_Ar.out"         , "-d 15 -c"),
        ("ref.Ar.no_sort.compress.out"              , "result.Ar.no_sort.compress.out"              , "Ar_Ar.out"         , "-d 15 --no-sort -c"),
        ("ref.N2.compress.out"                      , "result.N2.compress.out"                      , "N2_N2.out"         , "-d 15 -c"),
        ("ref.N2.compress.positronic.out"           , "result.N2.compress.positronic.out"           , "N2_N2.out"         , "-d 15 -pc"),
        ("ref.N2.compress.all.out"                  , "result.N2.compress.all.out"                  , "N2_N2.out"         , "-d 15 -ac"),
        ("ref.N2.no_sort.compress.out"              , "result.N2.no_sort.compress.out"              , "N2_N2.out"         , "-d 15 --no-sort -c"),
        ("ref.N2.no_sort.compress.positronic.out"   , "result.N2.no_sort.compress.positronic.out"   , "N2_N2.out"         , "-d 15 --no-sort -pc"),
        ("ref.N2.no_sort.compress.all.out"          , "result.N2.no_sort.compress.all.out"          , "N2_N2.out"         , "-d 15 --no-sort -ac"),
        ("ref.uo2.compress.out"                     , "result.uo2.compress.out"                     , "x2c_uo2_238.out"   , "-d 15 -c"),
        ("ref.uo2.no_sort.compress.out"             , "result.uo2.no_sort.compress.out"             , "x2c_uo2_238.out"   , "-d 15 --no-sort -c"),
        ("ref.ucl4.compress.out"                    , "result.ucl4.compress.out"                    , "x2c_ucl4.out"      , "-d 15 -c"),
        ("ref.ucl4.no_sort.compress.out"            , "result.ucl4.no_sort.compress.out"            , "x2c_ucl4.out"      , "-d 15 --no-sort -c"),
        ("ref.Cm3+_phen.compress.out"               , "result.Cm3+_phen.compress.out"               , "x2c_Cm3+_phen.out" , "-d 15 -c"),
    ]
    # fmt: on
)
def test_sum_dirac_dfcoeff_compress(ref_filename: str, result_filename: str, input_filename: str, options: str):
    test_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(test_path)
    print(test_path, " test start...")

    ref_filepath = os.path.join(test_path, "data", ref_filename)
    result_filepath = os.path.join(test_path, "results", result_filename)
    input_filepath = os.path.join(test_path, "data", input_filename)

    test_command = f"sum_dirac_dfcoef -i {input_filepath} -o {result_filepath} {options}"
    print(test_command)
    subprocess.run(
        test_command.split(),
        encoding="utf-8",
        check=True,
    )

    ref_file: List[List[str]] = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(ref_filepath).read().splitlines()))]
    out_file: List[List[str]] = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(result_filepath).read().splitlines()))]
    # File should have the same number of lines
    assert len(ref_file) == len(out_file), f"Number of lines in {ref_filename}(={len(ref_file)}) and {result_filename}(={len(out_file)}) are different."
    threshold: float = 1e-10
    checked = len(ref_file)
    for line_idx, (ref, out) in enumerate(zip(ref_file, out_file)):
        # 1st line has header information about eigenvalues
        # E1g closed 6 open 0 virtual 60 E1u closed 12 open 0 virtual 54
        if line_idx == 0:
            assert ref == out
            continue
        # ref[0]: irrep, ref[1]: energy order index in the irrep, ref[2]: energy, ref[3:]: Symmetry value and coefficient
        # (e.g.) E1u 19 -8.8824415703374 B3uUpx 49.999172476298732 B2uUpy 49.999172476298732
        assert ref[0] == out[0], f"irrep in line {line_idx} of {ref_filename} and {result_filename} are different."
        assert ref[1] == out[1], f"Energy order index in line {line_idx} of {ref_filename} and {result_filename} are different."
        assert abs(float(ref[2]) - float(out[2])) == pytest.approx(0, threshold), f"Energy in line {line_idx} of {ref_filename} and {result_filename} are different."
        for idx, (ref_val, out_val) in enumerate(zip(ref[3:], out[3:])):
            if idx % 2 == 0:
                assert ref_val == out_val, f"Symmetry value in line {line_idx} of {ref_filename} and {result_filename} are different."
            else:
                assert abs(float(ref_val) - float(out_val)) == pytest.approx(
                    0, threshold
                ), f"Contribution of the AO in the MO in line {line_idx} of {ref_filename} and {result_filename} are different."

    open(f"test.{input_filename}.log", "w").write(f"{checked} lines checked")


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
    ]
    # fmt: on
)
def test_sum_dirac_dfcoeff(ref_filename: str, result_filename: str, input_filename: str, options: str):
    test_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(test_path)
    print(test_path, " test start...")

    ref_filepath = os.path.join(test_path, "data", ref_filename)
    result_filepath = os.path.join(test_path, "results", result_filename)
    input_filepath = os.path.join(test_path, "data", input_filename)

    test_command = f"sum_dirac_dfcoef -i {input_filepath} -o {result_filepath} {options}"
    print(test_command)
    subprocess.run(
        test_command.split(),
        encoding="utf-8",
        check=True,
    )

    ref_file: List[List[str]] = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(ref_filepath).read().splitlines()))]
    out_file: List[List[str]] = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(result_filepath).read().splitlines()))]
    # File should have the same number of lines
    assert len(ref_file) == len(out_file), f"Number of lines in {ref_filename}(={len(ref_file)}) and {result_filename}(={len(out_file)}) are different."
    threshold: float = 1e-10
    checked = len(ref_file)
    for line_idx, (ref, out) in enumerate(zip(ref_file, out_file)):
        # 1st line has header information about eigenvalues
        # E1g closed 6 open 0 virtual 60 E1u closed 12 open 0 virtual 54
        if line_idx == 0:
            assert ref == out
            continue
        if len(ref) < 2 or len(out) < 2:
            checked -= 1
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
        assert ref_list_str == out_list_str, f"line {line_idx}: {ref_list_str} != {out_list_str}\nref: {ref_file[line_idx]}\nout:{out_file[line_idx]}"
        assert abs(ref_value - out_value) == pytest.approx(0, abs=threshold), f"line {line_idx}: {ref_value} != {out_value}\nref: {ref_file[line_idx]}\nout:{out_file[line_idx]}"
    open(f"test.{input_filename}.log", "w").write(f"{checked} lines checked")
