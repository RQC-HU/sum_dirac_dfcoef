import os
import pytest
import re
import subprocess
import sys


@pytest.mark.parametrize(
    "ref_filename, result_filename, input_filename, mol, options",
    # fmt: off
    [
        ("ref.Ar.out"                    , "result.Ar.out"                    , "Ar_Ar.out"       , "Ar"   , "-d 15"),
        ("ref.Ar.no_sort.out"            , "result.Ar.no_sort.out"            , "Ar_Ar.out"       , "Ar"   , "-d 15 --no-sort"),
        ("ref.Ar.compress.out"           , "result.Ar.compress.out"           , "Ar_Ar.out"       , "Ar"   , "-d 15 -c"),
        ("ref.Ar.no_sort.compress.out"   , "result.Ar.no_sort.compress.out"   , "Ar_Ar.out"       , "Ar"   , "-d 15 --no-sort -c"),
        ("ref.uo2.out"                   , "result.uo2.out"                   , "x2c_uo2_238.out" , "UO2"  , "-d 15"),
        ("ref.uo2.no_sort.out"           , "result.uo2.no_sort.out"           , "x2c_uo2_238.out" , "UO2"  , "-d 15 --no-sort"),
        ("ref.uo2.compress.out"          , "result.uo2.compress.out"          , "x2c_uo2_238.out" , "UO2"  , "-d 15 -c"),
        ("ref.uo2.no_sort.compress.out"  , "result.uo2.no_sort.compress.out"  , "x2c_uo2_238.out" , "UO2"  , "-d 15 --no-sort -c"),
        ("ref.ucl4.out"                  , "result.ucl4.out"                  , "x2c_ucl4.out"    , "UCl4" , "-d 15"),
        ("ref.ucl4.no_sort.out"          , "result.ucl4.no_sort.out"          , "x2c_ucl4.out"    , "UCl4" , "-d 15 --no-sort"),
        ("ref.ucl4.compress.out"         , "result.ucl4.compress.out"         , "x2c_ucl4.out"    , "UCl4" , "-d 15 -c"),
        ("ref.ucl4.no_sort.compress.out" , "result.ucl4.no_sort.compress.out" , "x2c_ucl4.out"    , "UCl4" , "-d 15 --no-sort -c")
    ]
    # fmt: on
)
def test_sum_dirac_dfcoeff(ref_filename: str, result_filename: str, input_filename: str, mol: str, options: str):

    script_name = "../sum_dirac_dfcoef"
    test_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(test_path)
    print(test_path, " test start...")

    ref_filepath = os.path.join(test_path, ref_filename)
    result_filepath = os.path.join(test_path, result_filename)
    script_filepath = os.path.join(test_path, script_name)

    test_command = f"{script_filepath} -m {mol} -i {input_filename} {options}"
    print(test_command)
    with open(result_filepath, "w") as file_output:
        process = subprocess.run(
            test_command,
            shell=True,
            stdout=file_output,
            encoding="utf-8",
            stderr=file_output,
        )
        if process.returncode != 0:
            sys.exit(f"{test_command} failed with return code {process.returncode}")

    ref_file: "list[list[str]]" = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(ref_filepath, "r").read().splitlines()))]
    out_file: "list[list[str]]" = [re.split(" +", line.rstrip("\n")) for line in list(filter(lambda val: val != "", open(result_filepath, "r").read().splitlines()))]
    # File should have the same number of lines
    assert len(ref_file) == len(out_file), f"Number of lines in {ref_filename}(={len(ref_file)}) and {result_filename}(={len(out_file)}) are different."
    threshold: float = 1e-10
    checked = len(ref_file)
    for line_idx, (ref, out) in enumerate(zip(ref_file, out_file)):
        if len(ref) < 2 or len(out) < 2:
            checked -= 1
            continue
        if "%" in ref[-1]:
            ref_value = float(ref[-2])
            out_value = float(out[-2])
        else:
            ref_value = float(ref[-1])
            out_value = float(out[-1])
        assert abs(ref_value - out_value) == pytest.approx(0, abs=threshold), f"line {line_idx}: {ref_value} != {out_value}\nref: {ref_file[line_idx]}\nout:{out_file[line_idx]}"
    open(f"test.{mol}.log", "w").write(f"{checked} lines checked")
