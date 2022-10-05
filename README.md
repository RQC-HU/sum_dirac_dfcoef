# SUMMARIZE DIRAC DFCOEF COEFFICIENTS

[![summarize_dirac_coefficients_test](https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/actions/workflows/ci.yml/badge.svg)](https://github.com/kohei-noda-qcrg/summarize_dirac_dfcoef_coefficients/actions/workflows/ci.yml)


This program summarizes the coefficients from a DIRAC output file using the *PRIVEC and .VECPRI options.  
(c.f. http://www.diracprogram.org/doc/master/manual/analyze/privec.html)


## Requirements

- [Python](https://python.org) (version ≧ 3.6)
  - If you don't know how to install python, I recommend you to use [pyenv](https://github.com/pyenv/pyenv)

## Usage

You can use this program with the following command!

```sh
/path/to/sum_dirac_dfcoef -f OUPUT_FILE_PATH -m MOLECULE_NAME
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

## Options

optional arguments (-f and -mol are required)

- -h, --help

  show this help message and exit  

- -f FILE, --file FILE

  (required) file name of DIRAC output

- -m MOL, --mol MOL

  (required) molecule specification. Write the molecular formula (e.g. Cu2O). *** DON'T write the rational formula (e.g. CH3OH) ***

- -t THRESHOLD, --threshold THRESHOLD

  threshold. Default: 0.1 % (e.g) --threshold=0.5 => print orbital with more than 0.5 % contribution
