# If user uses this package as a script, the following code is executed. Otherwise(used as a module), __init__.py is executed.
import sys

if __name__ == "__main__":
    from sum_dirac_dfcoef.sum_dirac_dfcoef import main

    sys.exit(main())
