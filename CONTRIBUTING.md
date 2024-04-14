## Contributing

First, thank you for considering contributing to this project!
You can contribute by opening a pull request or an issue.

- If you want to improve this project, please [fork this repository](https://github.com/RQC-HU/sum_dirac_dfcoef/fork) and [create a pull request](https://github.com/RQC-HU/sum_dirac_dfcoef/compare).
- If you find a bug, please [create an issue](https://github.com/RQC-HU/sum_dirac_dfcoef/issues/new).

## Setup development environment

  ```sh
  git clone https://github.com/yourname/sum_dirac_dfcoef.git
  cd sum_dirac_dfcoef
  pip install -e .[dev]
  ```

## CODING RULE

- This program writes header info to the output file if -c or -g option was set.
  - The format of the header information must follow the format below!
    - One blank line is required between header and data.
    - The first line must be key/value pairs only, values must not contain spaces.
    - If you want to add multiple values of header info, you must add them as a new line at the end of the header area.
    - If the value is to be used for [dcaspt2_input_generator](https://github.com/RQC-HU/dcaspt2_input_generator), please add the code to read the value to dcaspt2_input_generator as well
