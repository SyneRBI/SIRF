# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2018 - 2019 University College London
#
# This is software developed for the Collaborative Computational
# Project in Positron Emission Tomography and Magnetic Resonance imaging
# (http://www.ccppetmr.ac.uk/).
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# Imports
import os
import sys
import time
import numpy

import sirf.STIR as pet
import sirf.Gadgetron as mr
import sirf.Reg as reg
from sirf.Utilities import error

# Paths
SIRF_PATH = os.environ.get('SIRF_PATH')


def try_stirtonifti(nifti_filename):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting STIR to Nifti test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Load the image as a NiftiImageData3D
    image_nifti = reg.NiftiImageData3D(nifti_filename)

    # Read as STIRImageData, convert to NiftiImageData3D and save to file
    image_stir = pet.ImageData(nifti_filename)
    image_nifti_from_stir = reg.NiftiImageData3D(image_stir)
    image_nifti_from_stir.write('results/stir_to_nifti.nii',image_nifti.get_original_datatype())

    # Compare the two
    if image_nifti != image_nifti_from_stir:
        raise AssertionError("Conversion from STIR to Nifti failed.")

    # Resample and then check that voxel values match
    resample = reg.NiftyResample()
    resample.set_floating_image(image_stir) 
    resample.set_reference_image(image_nifti) 
    resample.set_interpolation_type_to_nearest_neighbour()
    resample.process()

    # as_array() of both original images should match
    if not numpy.array_equal(image_nifti.as_array(),resample.get_output().as_array()):
        raise AssertionError("as_array() of sirf.Reg.NiftiImageData and resampled sirf.STIR.ImageData are different.")

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished STIR to Nifti test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def try_gadgetrontonifti(nifti_filename, mr_recon_h5_filename):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Gadgetron to Nifti test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Read ISMRMRD image
    ismrmrd_im = mr.ImageData(mr_recon_h5_filename)

    # Convert ISMRMRD image to nifti
    nifti_from_ismrmrd = reg.ImageData(ismrmrd_im)

    # Read vendor-reconstructed image
    dicom_im = reg.ImageData(nifti_filename)

    # Standardise to remove scaling problems
    dicom_im.standardise()
    nifti_from_ismrmrd.standardise()

    # Compare the two. Since the images are being reconstructed independently, there is no
    # guarantee they will perfectly match. So we need an data-dependent acceptance threshold.
    if not reg.ImageData.are_equal_to_given_accuracy(dicom_im, nifti_from_ismrmrd, 165.):
        raise AssertionError("Conversion from ISMRMRD to Nifti failed")

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Gadgetron to Nifti test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# complex resample
def try_complex_resample(raw_mr_filename):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting complex resampling test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    raw_mr = mr.AcquisitionData(raw_mr_filename)

    recon_gadgets = ['NoiseAdjustGadget',
                     'AsymmetricEchoAdjustROGadget',
                     'RemoveROOversamplingGadget',
                     'AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)',
                     'BucketToBufferGadget(split_slices=true, verbose=false)',
                     'SimpleReconGadget',
                     'ImageArraySplitGadget']

    recon = mr.Reconstructor(recon_gadgets)
    recon.set_input(raw_mr)
    recon.process()

    ismrmrd_im = recon.get_output()

    if ismrmrd_im.is_real():
        raise AssertionError("Expected output of reconstruction to be complex")

    # Complex component may be empty, so use imag = real / 2
    image_data_arr = ismrmrd_im.as_array()
    image_data_arr.imag = image_data_arr.real / 2
    ismrmrd_im.fill(image_data_arr)

    # Convert the complex image to two niftis
    [real, imag] = reg.ImageData.construct_from_complex_image(ismrmrd_im)
    real.write("results/real")
    imag.write("results/imag")

    # Create affine transformation
    tm = reg.AffineTransformation()
    tm_ = tm.as_array()
    tm_[0][3] = 2.
    tm = reg.AffineTransformation(tm_)

    # Resample the complex data
    res_complex = reg.NiftyResample()
    res_complex.set_reference_image(ismrmrd_im)
    res_complex.set_floating_image(ismrmrd_im)
    res_complex.set_interpolation_type_to_linear()
    res_complex.add_transformation(tm)
    forward_cplx_sptr = res_complex.forward(ismrmrd_im)
    adjoint_cplx_sptr = res_complex.adjoint(ismrmrd_im)

    # Get the output
    [forward_cplx_real, forward_cplx_imag] = \
        reg.ImageData.construct_from_complex_image(forward_cplx_sptr)
    [adjoint_cplx_real, adjoint_cplx_imag] = \
        reg.ImageData.construct_from_complex_image(adjoint_cplx_sptr)

    forward_cplx_real.write("results/forward_cplx_real")
    forward_cplx_imag.write("results/forward_cplx_imag")
    adjoint_cplx_real.write("results/adjoint_cplx_real")
    adjoint_cplx_imag.write("results/adjoint_cplx_imag")

    # Now resample each of the components individually
    res_real = reg.NiftyResample()
    res_real.set_reference_image(real)
    res_real.set_floating_image(real)
    res_real.set_interpolation_type_to_linear()
    res_real.add_transformation(tm)
    forward_real = res_real.forward(real)
    adjoint_real = res_real.adjoint(real)

    res_imag = reg.NiftyResample()
    res_imag.set_reference_image(imag)
    res_imag.set_floating_image(imag)
    res_imag.set_interpolation_type_to_linear()
    res_imag.add_transformation(tm)
    forward_imag = res_imag.forward(imag)
    adjoint_imag = res_imag.adjoint(imag)

    forward_real.write("results/forward_real")
    forward_imag.write("results/forward_imag")
    adjoint_real.write("results/adjoint_real")
    adjoint_imag.write("results/adjoint_imag")

    # Compare that the real and imaginary parts match regardless
    # of whether they were resampled separately or together.
    if forward_real != forward_cplx_real or forward_imag != forward_cplx_imag:
        raise AssertionError("NiftyResample::forward failed for complex data")
    if adjoint_real != adjoint_cplx_real or adjoint_imag != adjoint_cplx_imag:
        raise AssertionError("NiftyResample::adjoint failed for complex data")

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished complex resampling test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
    raw_mr_filename = SIRF_PATH + "/data/examples/MR/grappa2_1rep.h5"
    if os.path.isfile(SIRF_PATH + "/data/examples/MR/zenodo/dicom_as_nifti.nii"):
        nifti_filename = SIRF_PATH + "/data/examples/MR/zenodo/dicom_as_nifti.nii"
        mr_recon_h5_filename = SIRF_PATH + "/data/examples/MR/zenodo/SIRF_recon.h5"
    else:
        nifti_filename = SIRF_PATH + "/data/examples/Registration/test2.nii.gz"

    try_stirtonifti(nifti_filename)
    if mr_recon_h5_filename:
        try_gadgetrontonifti(nifti_filename, mr_recon_h5_filename)
    try_complex_resample(raw_mr_filename)


if __name__ == "__main__":
    try:
        test()
    except:
        raise error("Error encountered.")
