'''Regsitration of SIRF images.

Usage:
  registration [--spm] [--help]

Options:
  --spm  test spm functionality
'''

# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2018 - 2020 University College London
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
import numpy as np
import sirf.Reg
from pUtilities import *

# Paths
SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration'
output_prefix = os.getcwd() + '/results/python_'

# Input filenames
ref_aladin_filename = examples_path + "/test.nii.gz"
spm_working_folder = output_prefix + "spm_working_folder"
spm_working_folder2 = output_prefix + "spm_working_folder2"
ref_aladin = sirf.Reg.NiftiImageData3D(ref_aladin_filename)
spm_to_register_ref = output_prefix + "spm_to_register_ref.nii"
spm_to_register_flo = output_prefix + "spm_to_register_flo.nii"


# SPM
def try_spm():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting SPM test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Resample an image with NiftyResample. Register SPM, check the result

    # TM
    translations = np.array([5.,  4., -5.], dtype=np.float32)
    euler_angles = np.array([5., -2., -3.], dtype=np.float32)

    tm = sirf.Reg.AffineTransformation(translations, euler_angles)

    niftyreg_resampler = sirf.Reg.NiftyResample()
    niftyreg_resampler.set_padding_value(0.)
    niftyreg_resampler.set_reference_image(ref_aladin)
    niftyreg_resampler.set_floating_image(ref_aladin)
    niftyreg_resampler.add_transformation(tm)
    niftyreg_resampler.set_interpolation_type_to_linear()
    floating = niftyreg_resampler.forward(ref_aladin)

    # Register with SPM
    spm_reg = sirf.Reg.SPMRegistration()
    spm_reg.set_reference_image(ref_aladin)
    spm_reg.add_floating_image(floating)
    spm_reg.add_floating_image(floating)
    spm_reg.set_working_folder(spm_working_folder)
    spm_reg.set_working_folder_file_overwrite(True)
    spm_reg.set_delete_temp_files(False)
    spm_reg.process()
    spm_tm = spm_reg.get_transformation_matrix_forward(1)
    spm_inv_tm = spm_tm.get_inverse()

    # Check tm roughly equals inverse TM of the resampler
    estimated_euler_angles = spm_inv_tm.get_Euler_angles()
    estimated_translations = spm_inv_tm.as_array()[:3, 3]

    input_euler_angles = tm.get_Euler_angles()
    input_translations = translations

    diff_euler_angles = 100. * (input_euler_angles - estimated_euler_angles) / input_euler_angles
    diff_translations = 100. * (input_translations - estimated_translations) / input_translations

    print("Input Euler angles:              " + np.array2string(input_euler_angles))
    print("Estimated Euler angles:          " + np.array2string(estimated_euler_angles))
    print("Percentage diff in Euler angles: " + np.array2string(diff_euler_angles))
    print("Input translations:              " + np.array2string(input_translations))
    print("Estimated translations:          " + np.array2string(estimated_translations))
    print("Percentage diff in translations: " + np.array2string(diff_translations))

    # Check differences are less than 1%
    if any(diff_euler_angles > 1):
        raise AssertionError("SPM registration failed (angles).")
    if any(diff_translations > 1):
        raise AssertionError("SPM registration failed (translations).")

    if spm_reg.get_output(1) != ref_aladin:
        raise AssertionError("SPM registration failed (image difference).")

    ref_aladin.write(spm_to_register_ref)
    floating.write(spm_to_register_flo)

    # Try to register via filename
    spm_reg2 = sirf.Reg.SPMRegistration()
    spm_reg2.set_reference_image_filename(spm_to_register_ref)
    spm_reg2.add_floating_image_filename(spm_to_register_flo)
    spm_reg2.add_floating_image_filename(spm_to_register_flo)
    spm_reg2.set_working_folder(spm_working_folder2)
    spm_reg2.set_working_folder_file_overwrite(True)
    spm_reg2.set_delete_temp_files(False)
    spm_reg2.process()

    for i in range(0, 2):
        spm_reg2.get_output(i).write(output_prefix + "spm_out_" + str(i))
        spm_reg2.get_displacement_field_forward(i).write(output_prefix + "spm_disp_fwd_" + str(i))
        spm_reg2.get_displacement_field_inverse(i).write(output_prefix + "spm_disp_inv_" + str(i))
        spm_reg2.get_deformation_field_forward(i).write(output_prefix + "spm_def_fwd_" + str(i))
        spm_reg2.get_deformation_field_inverse(i).write(output_prefix + "spm_def_inv_" + str(i))
        spm_reg2.get_transformation_matrix_forward(i).write(output_prefix + "spm_tm_fwd_" + str(i))
        spm_reg2.get_transformation_matrix_inverse(i).write(output_prefix + "spm_tm_inv_" + str(i))

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished SPM test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
    try_spm()


if __name__ == "__main__":
    try:
        test()
    except:
        raise error("Error encountered.")
