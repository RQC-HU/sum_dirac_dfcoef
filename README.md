# sum_dirac_dfcoef : SUMMARIZE DIRAC DFCOEF COEFFICIENTS

[![sum_dirac_dfcoef_test](https://github.com/RQC-HU/sum_dirac_dfcoef/actions/workflows/test.yml/badge.svg)](https://github.com/RQC-HU/sum_dirac_dfcoef/actions/workflows/test.yml)

This program provides a utility to summarize the contribution of each atomic orbital per kramers pair from the [DIRAC](http://diracprogram.org/doku.php) output file that the [.ANALYZE](https://diracprogram.org/doc/master/manual/dirac.html#analyze) and [.PRIVEC](https://diracprogram.org/doc/master/manual/analyze.html#privec) options are activated.

## Requirements

- [Python](https://python.org) (version ≧ 3.7)

## Install

```sh
pip install -U sum_dirac_dfcoef
```

## Usage

### Linux, macOS

You can use this program with the following command!

```sh
# Output to sum_dirac_dfcoef.out
sum_dirac_dfcoef -i DIRAC_OUPUT_FILE_PATH
# Specify output file name with -o option
sum_dirac_dfcoef -i DIRAC_OUPUT_FILE_PATH -o OUTPUT_FILE_NAME
```

(e.g.)

```sh
sum_dirac_dfcoef -i x2c_uo2_238.out
```

### Windows

If you want to use this program on Windows, you can use it with the following command.

```sh
sum_dirac_dfcoef.exe -i DIRAC_OUPUT_FILE_PATH
# or
python -m sum_dirac_dfcoef -i DIRAC_OUPUT_FILE_PATH
```

A part of x2c_uo2_238.out (DIRAC output file, ... represents an omission)

```out
...
    **************************************************************************
    ****************************** Vector print ******************************
    **************************************************************************



   Coefficients from DFCOEF
   ------------------------



                                Fermion ircop E1g
                                -----------------


* Electronic eigenvalue no. 17: -5.1175267254674
====================================================
       1  L Ag U  s      -0.0000003723  0.0000000000  0.0000000000  0.0000000000
       2  L Ag U  s      -0.0000008538  0.0000000000  0.0000000000  0.0000000000
       3  L Ag U  s      -0.0000014888  0.0000000000  0.0000000000  0.0000000000
       4  L Ag U  s      -0.0000025924  0.0000000000  0.0000000000  0.0000000000
       5  L Ag U  s      -0.0000043736  0.0000000000  0.0000000000  0.0000000000
       6  L Ag U  s      -0.0000074960  0.0000000000  0.0000000000  0.0000000000
...

*****************************************************
********** E N D   of   D I R A C  output  **********
*****************************************************
...
```

A part of the result (... represents an omission)

```out
NO_HEADERINFO: This output cannot be used for the dcaspt2_input_generator program.
NO_HEADERINFO: This output cannot be used for the dcaspt2_input_generator program.

Electronic no. 19 E1u -8.88244
B3uUpx        49.99917 %
B2uUpy        49.99917 %

Electronic no. 20 E1u -8.86075
B1uUpz        66.76661 %
B3uUpx        16.05235 %
B2uUpy        16.05235 %
B1uOs(1)       0.54741 %
B1uOs(2)       0.54741 %

Electronic no. 17 E1g -5.11753
B2gUdxz       35.98781 %
B3gUdyz       35.98781 %
AgUdzz        18.54868 %
AgUdxx         4.63717 %
AgUdyy         4.63717 %
AgUs           0.13729 %
...
```

If you use -c or --compress option, you can get a compressed result like this.(one line per kramers pair)

```out
electron_num 106 point_group D2h moltra_scheme default
E1g 16..85 E1u 11..91 
E1g closed 52 open 0 virtual 268 E1u closed 54 open 0 virtual 314 

E1u 19 -8.88244 B3uUpx 49.99917 B2uUpy 49.99917
E1u 20 -8.86075 B1uUpz 66.76661 B3uUpx 16.05235 B2uUpy 16.05235 B1uOs(1) 0.54741 B1uOs(2) 0.54741
E1g 17 -5.11753 B2gUdxz 35.98781 B3gUdyz 35.98781 AgUdzz 18.54868 AgUdxx 4.63717 AgUdyy 4.63717 AgUs 0.13729
...
```

This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.

## dcaspt2_input_generator

- If you want to use the result of this program as input to [dcaspt2_input_generator](https://github.com/RQC-HU/dcaspt2_input_generator), we recommend you to use -g or --for-generator option.

## Options

optional arguments (--input is required)

- -h, --help

  show this help message and exit

- -i "INPUT", --input "INPUT"

  (required) file path of DIRAC output.  
  Please quote if the path include spaces.

- -o "OUTPUT", --output "OUTPUT"

  File path of sum_dirac_dfcoef output.  
  Default: sum_dirac_dfcoef.out.  
  Please quote if the path include spaces.  

- -g, --for-generator

  Automatically set the arguments for dcaspt2_input_generator.  
  This option is useful when you want to use the result of this program as input to dcaspt2_input_generator.  
  This option is equivalent to set -c/--compress and not set -p/--positronic and --no-scf options.

- -j [PARALLEL], --parallel [PARALLEL]

  Number of parallel processes.  
  Default: 1 (single process).  
  If you set -j option without argument, the number of parallel processes is set to the number of CPU cores(=os.cpu_count()).

- -c, --compress

  Compress output.  
  Display coefficients on one line for each kramers pair.  
  This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.

- -t THRESHOLD, --threshold THRESHOLD

  threshold. Default: 0.1 %  
  (e.g) --threshold=0.1 → only print atomic orbitals with more than 0.1 % contribution

- -d DECIMAL, --decimal DECIMAL

  Set the decimal places.  
  Default: 5  
  (e.g) --decimal=3 → print orbital with 3 decimal places (0.123, 2.456, ...).  
  range: 1-15

- --ignore-atom-num

  Ignore the atom number label.  
  This option is useful when you want to sum the coefficients of the same atomic orbital except for the atom number label.

- --ignore-ml

  Ignore the magnetic quantum number label.  
  This option is useful when you want to sum the coefficients of the same atomic orbital except for the magnetic quantum number label.

- --ignore-sym

  Ignore symmetry label (e.g. Ag, B3g, ...).  
  This option is useful when you want to sum the coefficients of the same atomic orbital except for the symmetry label.

- -a, --all-write

  Print all kramers pairs(Positronic and Electronic).

- -p, --positronic-write

  Print only Positronic kramers pairs.  
  The output with this option cannot be used as input to dcaspt2_input_generator.

- -v, --version

  Print version and exit

- --no-scf

  If you don't activate .SCF keyword in your DIRAC input file, you must use this option.  
  But you cannot use the output using this option to dcaspt2_input_generator program.

- --debug

  print debug output (Normalization constant, Sum of kramers pair coefficient ...)

- --no-sort

  Don't sort the output by kramers pair energy

## Development

- Thank you for considering contributing to this project!
- Please read [CONTRIBUTING.md](./CONTRIBUTING.md) before you start contributing.

## LICENSE

- LGPL v2.1

## Maintainer

- Minori Abe

## Author

- Kohei Noda
