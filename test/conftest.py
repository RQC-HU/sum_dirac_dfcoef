import os
import subprocess
from pathlib import Path


def pytest_sessionstart():
    cur_dir = Path.cwd()
    # Change the current directory to the root directory of the package
    test_dir = Path(__file__).resolve().parent
    results_dir = Path.joinpath(test_dir, "results")
    if results_dir.exists():
        # Remove *.out files in the results directory
        for file in results_dir.glob("*.out"):
            os.remove(file)
    os.chdir(Path.joinpath(test_dir, ".."))
    cmd = "python3 -m pip install -e .[test]"
    subprocess.run(cmd.split(), check=True)
    # Change the current directory to the original directory
    os.chdir(cur_dir)
