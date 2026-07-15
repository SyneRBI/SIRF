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

### Examples data path

The resolution order for the examples data path is:

1. The directory pointed by `${SIRF_DATA_PATH}/examples` if the environment variable `SIRF_DATA_PATH` is set by the user.
2. The directory pointed by the CMake variable `SIRF_EXAMPLES_DATA_PATH` if it is set by the user during building SIRF.
3. The directory `${CMAKE_INSTALL_PREFIX}/share/SIRF-<version_major>.<version_minor>`
The user may need options 1 or 2 e.g. if the examples input data is in a write-protected location, and the examples are designed to send output there, hence a copy to writable location is needed.
