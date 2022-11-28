# SUMMARIZE DIRAC DFCOEF COEFFICIENTS

[![summarize_dirac_coefficients_test](https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/actions/workflows/ci.yml/badge.svg)](https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/actions/workflows/ci.yml)


This program summarizes the coefficients from a DIRAC output file using the *PRIVEC and .VECPRI options.  
(c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)


## Requirements

- [Python](https://python.org) (version ≧ 3.6)
  - If you don't know how to install python, I recommend you to use [pyenv](https://github.com/pyenv/pyenv)

## Download

```sh
git clone https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients.git
```

## Usage

You can use this program with the following command!

```sh
/path/to/sum_dirac_dfcoef -f DIRAC_OUPUT_FILE_PATH -m MOLECULE_NAME
```

(e.g.)

```sh
./sum_dirac_dfcoef -f x2c_uo2_238.out -m UO2
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
Electronic no. 17 E1g -5.1175267254674
AgUs       :   0.13728512177623234 %
AgUdxx     :   4.63716965312627316 %
AgUdyy     :   4.63716965312627316 %
AgUdzz     :  18.54867861134004769 %
B2gUdxz    :  35.98780872755493476 %
B3gUdyz    :  35.98780872755493476 %
Normalization constant is 0.5059238161886942
sum of coefficient 0.9999999999999996

Electronic no. 18 E1g -5.1107907830158
AgUdxx     :  12.47047856323067627 %
AgUdyy     :  12.47047856323067627 %
B1gUdxy    :  49.88191426953538610 %
B2gUdxz    :  12.58821772533592309 %
B3gUdyz    :  12.58821772533592309 %
Normalization constant is 0.5177168489972134
...
```

If you use -c or --compress option, you can get a compressed result like this.(one line per MO)

This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.

```out
E1g 17 -5.1175267254674 AgUs 0.137285121776232 AgUdxx 4.637169653126273 AgUdyy 4.637169653126273 AgUdzz 18.548678611340048 B2gUdxz 35.987808727554935 B3gUdyz 35.987808727554935
E1g 18 -5.1107907830158 AgUdxx 12.470478563230676 AgUdyy 12.470478563230676 B1gUdxy 49.881914269535386 B2gUdxz 12.588217725335923 B3gUdyz 12.588217725335923
E1g 19 -4.8038359333701 AgUs 0.134289378478735 AgUdxx 6.527712047428291 AgUdyy 6.527712047428291 AgUdzz 26.110848178313287 B2gUdxz 30.337915669150643 B3gUdyz 30.337915669150643
...
```


## Options

optional arguments (-f and -mol are required)

- -h, --help

  show this help message and exit  

- -f FILE, --file FILE

  (required) file name of DIRAC output

- -m MOL, --mol MOL

  (required) molecule specification. Write the molecular formula (e.g. Cu2O). *** DON'T write the rational formula (e.g. CH3OH) ***

- -c, --compress
  Compress output. Display all coefficients on one line for each MO.
  This options is useful when you want to use the result in a spreadsheet like Microsoft Excel.

- -t THRESHOLD, --threshold THRESHOLD

  threshold. Default: 0.1 % (e.g) --threshold=0.1 => print orbital with more than 0.1 % contribution

- -d DECIMAL, --decimal DECIMAL
Set the decimal places. Default: 5 (e.g) --decimal=3 => print orbital with 3 decimal places (0.123, 2.456, ...). range: 1-15
- --debug               print debug output (Normalization constant, Sum of MO coefficient)