#!/usr/bin/env python3

import shutil
from pathlib import Path


def main():
    """This script is used to update the ref.*.out files in the test/references directory
    with the result.*.out files in the test/results directory.

    This script is used when the results of the program are changed.
    """
    this_file_dir = Path(__file__).resolve().parent
    results_dir = this_file_dir / "results"
    references_dir = this_file_dir / "references"
    result_list = list(results_dir.glob("result.*.out"))
    for result in result_list:
        # replace only the first occurrence of "result" with "reference"
        ref_name = result.name.replace("result", "ref", 1)
        ref_path = references_dir / ref_name
        shutil.copy2(result, ref_path)


if __name__ == "__main__":
    main()
