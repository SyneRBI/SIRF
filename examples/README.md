This directory contains examples and demos on how to use SIRF. 

Please also our Jupyter notebooks in the [SIRF Exercises](https://github.com/SyneRBI/SIRF-Exercises/#readme)
which are more extensive and allow you to learn step-by-step.


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

The contents of subfolder named Python are described [here](Python/README.md). The contents of subfolder named Matlab are of similar nature but are not documented as we no longer support Matlab.

## Environment variables

The correct setup of SIRF requires defining the environment variables `SIRF_PATH` and `SIRF_INSTALL_PATH`.

If you built SIRF with the SIRF-Superbuild, these environment variables can be sourced from the `env_sirf.*` files, e.g.
on Unix/Linux [`env_sirf.sh`](https://github.com/SyneRBI/SIRF-SuperBuild/blob/master/env_sirf.sh.in).

Alternatively,  `SIRF_PATH` should be set by the user to the full path to the root folder containing the SIRF source,
and `SIRF_INSTALL_PATH` to the full path to the folder where the user installed SIRF binaries and other files (i.e.
set by `CMAKE_INSTALL_PREFIX`).

Additionally, the user may define the variable `SIRF_DATA_PATH` to point to an alternative location. For instance, when the
installation files are in a write-protected location, the user can make a copy elsewhere and point there.

<!--
(NOTE: even on Windows, you must use `/` in paths, not `\\`.)
-->

### Examples data path

The resolution order for the examples data path is:

1. The directory pointed by `${SIRF_DATA_PATH}` if it is set by the user.
2. The directory `${SIRF_INSTALL_PATH}/share/SIRF-<version_major>.<version_minor>` if `SIRF_INSTALL_PATH` is set.
3. The directory `${SIRF_PATH}/data/examples` if `SIRF_PATH` is set.
