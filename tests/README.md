# Tests in this folder

## NiftyPET_test.sh

This bash script compares the STIR-wrapped NiftyPET functionality to NiftyPET and, where possible, STIR.

### Usage

`NiftyPET_test.sh RAW_DATA_FOLDER TSTART TSTOP`

### Test content

From raw dicom data (`RAW_DATA_FOLDER`) and a given time window (`TSTART` and `TSTOP`), prompt sinograms are extracted with NiftyPET, STIR-NiftyPET and STIR. Randoms are estimated, and for NiftyPET and STIR-NiftyPET, a norm sinogram is extracted.

Forward and back projections are performed with NiftyPET and STIR-NiftyPET, first with no corrections, then with norm, randoms and attenuation taken into account.

Lastly, a single MLEM iteration is performed with NiftyPET and STIR-NiftyPET.

At each step, the results are compared to ensure consistency.

At the time of writing, NiftyPET requires `python2`, so two different python executables should be set as environmental variables - `NIFTYPET_PYTHON_EXECUTABLE` and `SIRF_PYTHON_EXECUTABLE`.

### Data

An example of test data can be found here - [https://doi.org/10.5281/zenodo.1472951](https://doi.org/10.5281/zenodo.1472951).

Since the data needs to be read by NiftyPET, it should be in the raw form of `.dcm`/`.bf`, and not `.l.hdr`/`.l`.