import os
import pytest
import re
import subprocess
import sys


def run_script_and_check(ref_filename: str, result_filename: str, input_filename: str, mol: str):

    script_name = "../sum_dirac_dfcoef.py"
    test_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(test_path)
    print(test_path, " test start...")

    ref_filepath = os.path.join(test_path, ref_filename)
    result_filepath = os.path.join(test_path, result_filename)
    script_filepath = os.path.join(test_path, script_name)

    test_command = f"python {script_filepath} -mol {mol} -f {input_filename}"
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
    open(f'test.{mol}.log', "w").write(f"{checked} lines checked")


def test_ucl4():
    ref_filename = "ref.ucl4.out"
    result_filename = "result.ucl4.out"
    input_filename = "x2c_ucl4.out"
    mol = "UCl4"
    run_script_and_check(ref_filename, result_filename, input_filename, mol)


def test_uo2():
    ref_filename = "ref.uo2.out"
    result_filename = "result.uo2.out"
    input_filename = "x2c_uo2_238.out"
    mol = "UO2"
    run_script_and_check(ref_filename, result_filename, input_filename, mol)


if __name__ == "__main__":
    test_uo2()
