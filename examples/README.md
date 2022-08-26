This directory contains examples and demos on how to use SIRF. 
<!--
Please also see the [related Wiki page](https://github.com/CCPPETMR/SIRF/wiki/Examples).
-->

Most of the demos use simplistic data for illustration. The images might therefore be
not very interesting either. The idea here is to illustrate what you can do with SIRF,
not to give you a polished script that you would want to use on a daily basis.

The demos display usually some images. Generally, these might not make too much sense
unless you have a look at the code.

The structure of the examples subfolder is as follows:

    examples
        Matlab
            MR
                Gadgetron
            PET
            PETMR
            Registration
        Python
            MR
                Gadgetron
            PET
            PETMR
            Registration
            SPECT
                interactive

The contents of subfolder named Python are described [here](https://github.com/SyneRBI/SIRF/blob/master/examples/Python/README.md). The contents of subfolder named Matlab are of similar nature but are not documented as we no longer support Matlab.

## Environment variables

The correct setup of SIRF requires defining the environment variables `SIRF_PATH` and `SIRF_INSTALL_PATH`.

On Unix/Linux, these environment variables are sourced from the [`env_sirf.sh`](https://github.com/SyneRBI/SIRF-SuperBuild/blob/master/env_sirf.sh.in) file.

On Windows, the value of `SIRF_PATH` should be set by the user to the full path to the root folder containing SIRF source, and the value of `SIRF_INSTALL_PATH` to the full path to the folder where the user wants SIRF binaries and other files created by the build installed.

Additionally, the user may define the variable `SIRF_DATA_PATH` whose value is the path to the data to be processed by SIRF.

### Examples data path

The resolution order for the examples data path is:

1. The directory pointed by `${SIRF_DATA_PATH}` if it is set by the user.
2. The directory `${SIRF_INSTALL_PATH}/share/SIRF-<version_major>.<version_minor>`.
