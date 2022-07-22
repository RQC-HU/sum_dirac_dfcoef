# SUMMARIZE DIRAC DFCOEF COEFFICIENTS

This program summarize the coefficients from DIRAC output file that *PRIVEC option was used.  
(c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)

## Requirements

- [Python](python.org) (version ≧ 3.6)
  - If you don't know how to install python, I recommend you to use [pyenv](https://github.com/pyenv/pyenv)

## Usage

You can use this program with the following command!

```sh
python summarize_dirac_dfcoef_coefficients.py -f OUPUT_FILE_PATH -mol MOLECULE_NAME
```

(e.g.)

```sh
python summarize_dirac_dfcoef_coefficients.py -f uo2_uo2.out -mol UO2
```

A part of uo2_uo2.out (DIRAC output file, ... represents an omission)

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
Electronic no. 17 E1g -5.1175267254674
AgUs       : 0.13728512177623234        %
AgUdxx     : 4.637169653126273          %
AgUdyy     : 4.637169653126273          %
AgUdzz     : 18.548678611340048         %
B2gUdxz    : 35.987808727554935         %
B3gUdyz    : 35.987808727554935         %
Normalization constant is 0.5059238161886942
sum of coefficient 0.9999999999999996

Electronic no. 18 E1g -5.1107907830158
AgUdxx     : 12.470478563230676         %
AgUdyy     : 12.470478563230676         %
B1gUdxy    : 49.881914269535386         %
B2gUdxz    : 12.588217725335923         %
B3gUdyz    : 12.588217725335923         %
Normalization constant is 0.5177168489972134
...
```

## Options

optional arguments (-f and -mol is required)

- -h, --help

  show this help message and exit  

- -f FILE, --file FILE

  (required) file name of DIRAC output

- -mol MOLECULE, --molecule MOLECULE

  (required) molecule specification. Write the molecular formula (e.g. Cu2O)

- -t THRESHOLD, --threshold THRESHOLD

  threshold. Default: 0.1 % (e.g) --threshold=0.5 => print orbital with more than 0.5 % contribution
