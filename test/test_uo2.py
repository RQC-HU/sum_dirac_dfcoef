import difflib
import os
import subprocess
import sys


def test_uo2():
    ref_filename = "ref.out"
    result_filename = "result.out"
    script_name = "../sum_dirac_dfcoef.py"

    test_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(test_path)
    print(test_path, " test start...")

    ref_filepath = os.path.join(test_path, ref_filename)
    result_filepath = os.path.join(test_path, result_filename)
    script_filepath = os.path.join(test_path, script_name)

    test_command = f"python {script_filepath} -mol UO2 -f x2c_uo2_238.out"
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

    ref_file = open(ref_filepath, "r")
    out_file = open(result_filepath, "r")

    diff = difflib.unified_diff(ref_file.readlines(), out_file.readlines())
    diff = "".join(diff)
    assert diff == ""


if __name__ == "__main__":
    test_uo2()
