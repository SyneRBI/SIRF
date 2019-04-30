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

import numpy as np
import pReg
from pUtilities import *

# Paths
SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration'
output_prefix = os.getcwd() + '/results/python_'

# Input filenames
ref_aladin_filename = examples_path + "/test.nii.gz"
flo_aladin_filename = examples_path + "/test2.nii.gz"
ref_f3d_filename = examples_path + "/mouseFixed.nii.gz"
flo_f3d_filename = examples_path + "/mouseMoving.nii.gz"
parameter_file_aladin = examples_path + "/paramFiles/niftyreg_aladin.par"
parameter_file_f3d = examples_path + "/paramFiles/niftyreg_f3d.par"

# Output filenames
save_nifti_image = output_prefix + "save_NiftiImageData.nii"
save_nifti_image_3d = output_prefix + "save_NiftiImageData3D.nii"
save_nifti_image_3d_tensor_not_split = output_prefix + "save_NiftiImageData3DTensor_not_split.nii"
save_nifti_image_3d_tensor_split = output_prefix + "save_NiftiImageData3DTensor_split_%s.nii"
save_nifti_image_3d_deformation_not_split = output_prefix + "save_NiftiImageData3DDeformation_not_split.nii"
save_nifti_image_3d_deformation_split = output_prefix + "save_NiftiImageData3DDeformation_split_%s.nii"
save_nifti_image_3d_displacement_not_split = output_prefix + "save_NiftiImageData3DDisplacement_not_split.nii"
save_nifti_image_3d_displacement_split = output_prefix + "save_NiftiImageData3DDisplacement_split_%s.nii"
save_nifti_image_upsample = output_prefix + "save_NiftiImageData_upsample.nii";
save_nifti_image_downsample = output_prefix + "save_NiftiImageData_downsample.nii";
save_nifti_image_up_downsample = output_prefix + "save_NiftiImageData_upsample_downsample.nii";
aladin_warped = output_prefix + "aladin_warped.nii"
f3d_warped = output_prefix + "f3d_warped.nii"
TM_forward = output_prefix + "TM_forward.txt"
TM_inverse = output_prefix + "TM_inverse.txt"
aladin_def_forward = output_prefix + "aladin_def_forward.nii"
aladin_def_inverse = output_prefix + "aladin_def_inverse_%s.nii"
aladin_disp_forward = output_prefix + "aladin_disp_forward.nii"
aladin_disp_inverse = output_prefix + "aladin_disp_inverse_%s.nii"
f3d_def_forward = output_prefix + "f3d_disp_forward.nii"
f3d_def_inverse = output_prefix + "f3d_disp_inverse_%s.nii"
f3d_disp_forward = output_prefix + "f3d_disp_forward.nii"
f3d_disp_inverse = output_prefix + "f3d_disp_inverse_%s.nii"

rigid_resample = output_prefix + "rigid_resample.nii"
nonrigid_resample_disp = output_prefix + "nonrigid_resample_disp.nii"
nonrigid_resample_def = output_prefix + "nonrigid_resample_def.nii"
output_weighted_mean = output_prefix + "weighted_mean.nii"
output_weighted_mean_def = output_prefix + "weighted_mean_def.nii"
output_float = output_prefix + "reg_aladin_float.nii"

ref_aladin = pReg.NiftiImageData3D(ref_aladin_filename)
flo_aladin = pReg.NiftiImageData3D(flo_aladin_filename)
ref_f3d = pReg.NiftiImageData3D(ref_f3d_filename)
flo_f3d = pReg.NiftiImageData3D(flo_f3d_filename)

# NiftiImageData
def try_niftiimage():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = pReg.NiftiImageData()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = pReg.NiftiImageData(ref_aladin_filename)

    # Save to file
    b.write(save_nifti_image)

    # Fill
    b.fill(100)

    # Get max
    if b.get_max() != 100:
        raise AssertionError('NiftiImageData fill()/get_max() failed.')

    # Get min
    if b.get_min() != 100:
        raise AssertionError('NiftiImageData fill()/get_min() failed.')

    # Deep copy
    d = b.deep_copy()
    if d.handle == b.handle:
        raise AssertionError('NiftiImageData deep_copy failed.')
    if d != b:
        raise AssertionError("NiftiImageData deep_copy failed.")

    # Addition
    e = d + d
    if abs(e.get_max() - 2 * d.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData __add__/get_max() failed.')

    # Subtraction
    e = d - d
    if abs(e.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData __sub__ failed.')

    # Sum
    if abs(e.get_sum()) > 0.0001:
        raise AssertionError('NiftiImageData get_sum() failed.')

    # Add num to image
    q = e + 1
    if q.get_max() != e.get_max() + 1:
        raise AssertionError('NiftiImageData __add__ val failed.')

    # Subtract num from image
    r = e - 1
    if r.get_max() != e.get_max() - 1:
        raise AssertionError('NiftiImageData __sub__ val failed.')

    # Multiply image by num
    s = e * 10
    if s.get_max() != e.get_max() * 10:
        raise AssertionError('NiftiImageData __mul__ val failed.')

    # Dimensions
    f = e.get_dimensions()
    if not np.array_equal(f, [3, 64, 64, 64, 1, 1, 1, 1]):
        raise AssertionError('NiftiImageData get_dimensions() failed.')

    # Get as array
    arr = d.as_array()
    if arr.max() != 100:
        raise AssertionError('NiftiImageData as_array().max() failed.')
    if arr.ndim != 3:
        raise AssertionError('NiftiImageData as_array() ndims failed.')
    if arr.shape != (64, 64, 64):
        raise AssertionError('NiftiImageData as_array().shape failed.')

    # Test saving to datatype
    ref_aladin.write(output_float, 16) # float
    ref_aladin_float = pReg.NiftiImageData3D(output_float)
    arr1 = ref_aladin.as_array()
    arr2 = ref_aladin_float.as_array()
    if not np.array_equal(arr1,arr2):
        raise AssertionError("NiftiImageData::write()/change_datatype() failed.")

    # Test print methods
    q.print_header()
    pReg.NiftiImageData.print_headers([q, s])

    # Crop image
    min_ = []
    max_ = []
    for i in range(0, 7):
        min_.append(0)
        max_.append(f[i+1] - 1)
    max_[2] = 62
    s = e
    s.crop(min_, max_)
    if s.as_array().shape != (64, 64, 63):
        raise AssertionError("NiftiImageData crop() failed.")

    # Get voxel sizes
    s = b.get_voxel_sizes()
    if not all(numpy.equal(s,numpy.array([0, 4.0625, 4.0625, 4.0625, 0, 0, 0, 0]))):
        raise AssertionError("NiftiImageData get_voxel_sizes() failed.")

    # Check upsampling/downsampling
    u = pReg.NiftiImageData(ref_aladin_filename);
    original_spacing    = u.get_voxel_sizes();
    original_spacing    = original_spacing[1:4];
    upsampled_spacing   = [original_spacing[0]/2, original_spacing[1]/4, original_spacing[2]];
    downsampled_spacing = [original_spacing[0]*2, original_spacing[1]*4, original_spacing[2]];
    # Downsample
    v = u.deep_copy();
    v.set_voxel_spacing(downsampled_spacing,3);
    v.write(save_nifti_image_downsample);
    # Upsample then downsample, check nothing has changed
    w = u.deep_copy();
    w.set_voxel_spacing(upsampled_spacing,0);
    w.write(save_nifti_image_upsample);
    x = w.deep_copy();
    x.set_voxel_spacing(original_spacing,0);
    x.write(save_nifti_image_up_downsample);
    pReg.NiftiImageData.print_headers([u, v, w, x]);
    if x != u:
        raise AssertionError('NiftiImageData::upsample()/downsample() failed.')


    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImageData test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImageData3D
def try_niftiimage3d():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData3D test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = pReg.NiftiImageData3D()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = pReg.NiftiImageData3D(ref_aladin_filename)

    # Save to file
    b.write(save_nifti_image_3d)

    # Fill
    b.fill(100)

    # Get max
    if b.get_max() != 100:
        raise AssertionError('NiftiImageData3D fill()/get_max() failed.')

    # Get min
    if b.get_min() != 100:
        raise AssertionError('NiftiImageData3D fill()/get_min() failed.')

    # Deep copy
    d = b.deep_copy()
    if d.handle == b.handle:
        raise AssertionError('NiftiImageData3D deep_copy failed.')
    if d != b:
        raise AssertionError("NiftiImageData3D deep_copy failed.")

    # Addition
    e = d + d
    if abs(e.get_max() - 2 * d.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3D __add__/get_max() failed.')

    # Subtraction
    e = d - d
    if abs(e.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3D __sub__ failed.')

    # Sum
    if abs(e.get_sum()) > 0.0001:
        raise AssertionError('NiftiImageData3D get_sum() failed.')

    # Dimensions
    f = e.get_dimensions()
    if not np.array_equal(f, [3, 64, 64, 64, 1, 1, 1, 1]):
        raise AssertionError('NiftiImageData3D get_dimensions() failed.')

    # Get as array
    arr = d.as_array()
    if arr.max() != 100:
        raise AssertionError('NiftiImageData3D as_array().max() failed.')
    if arr.ndim != 3:
        raise AssertionError('NiftiImageData3D as_array() ndims failed.')
    if arr.shape != (64, 64, 64):
        raise AssertionError('NiftiImageData3D as_array().shape failed.')

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImageData3D test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImageData3DTensor
def try_niftiimage3dtensor():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData3DTensor test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImageData3DTensor from NiftiImageData3D
    b = pReg.NiftiImageData3DTensor()
    b.create_from_3D_image(ref_aladin)

    # # Save to file
    b.write(save_nifti_image_3d_tensor_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_tensor_split)

    # Constructor from file
    c = pReg.NiftiImageData3DTensor(save_nifti_image_3d_tensor_not_split)

    # Fill
    c.fill(100)

    # Get max
    if c.get_max() != 100:
        raise AssertionError('NiftiImageData3DTensor fill()/get_max() failed.')

    # Get min
    if c.get_min() != 100:
        raise AssertionError('NiftiImageData3DTensor fill()/get_min() failed.')

    # Deep copy
    d = c.deep_copy()
    if d.handle == c.handle:
        raise AssertionError('NiftiImageData3DTensor deep_copy failed.')
    if d != c:
        raise AssertionError("NiftiImageData3DTensor deep_copy failed.")

    # Addition
    e = d + d
    if abs(e.get_max() - 2 * d.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DTensor __add__/get_max() failed.')

    # Subtraction
    e = d - d
    if abs(e.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DTensor __sub__ failed.')

    # Sum
    if abs(e.get_sum()) > 0.0001:
        raise AssertionError('NiftiImageData3DTensor get_sum() failed.')

    # Dimensions
    f = e.get_dimensions()
    if not np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]):
        raise AssertionError('NiftiImageData3DTensor get_dimensions() failed.')

    # Get as array
    arr = d.as_array()
    if arr.max() != 100:
        raise AssertionError('NiftiImageData3DTensor as_array().max() failed.')
    if arr.ndim != 5:
        raise AssertionError('NiftiImageData3DTensor as_array() ndims failed.')
    if arr.shape != (64, 64, 64, 1, 3):
        raise AssertionError('NiftiImageData3DTensor as_array().shape failed.')

    # Constructor from single components
    im1 = ref_aladin.deep_copy()
    im2 = ref_aladin.deep_copy()
    im3 = ref_aladin.deep_copy()
    im1.fill(30)
    im2.fill(20)
    im3.fill(-10)
    h = pReg.NiftiImageData3DTensor(im1, im2, im3)

    # Test flip components
    h.flip_component(0)
    if h.get_max() != 20:
        raise AssertionError("NiftiImageData3DTensor flip_component() failed.")
    if h.get_min() != -30:
        raise AssertionError("NiftiImageData3DTensor flip_component() failed.")


    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImageData3DTensor test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImageData3DDisplacement
def try_niftiimage3ddisplacement():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData3DDisplacement test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImageData3DDisplacement from NiftiImageData3D
    b = pReg.NiftiImageData3DDisplacement()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.write(save_nifti_image_3d_displacement_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_displacement_split)

    # Constructor from file
    c = pReg.NiftiImageData3DDisplacement(save_nifti_image_3d_displacement_not_split)

    # Constructor from 3x3D
    d = pReg.NiftiImageData3DDisplacement(ref_aladin, ref_aladin, ref_aladin)

    # Fill
    c.fill(100)

    # Get max
    if c.get_max() != 100:
        raise AssertionError('NiftiImageData3DDisplacement fill()/get_max() failed.')

    # Get min
    if c.get_min() != 100:
        raise AssertionError('NiftiImageData3DDisplacement fill()/get_min() failed.')

    # Deep copy
    d = c.deep_copy()
    if d.handle == c.handle:
        raise AssertionError('NiftiImageData3DDisplacement deep_copy failed.')
    if d != c:
        raise AssertionError("NiftiImageData3DDisplacement deep_copy failed.")

    # Addition
    e = d + d
    if abs(e.get_max() - 2 * d.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DDisplacement __add__/get_max() failed.')

    # Subtraction
    e = d - d
    if abs(e.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DDisplacement __sub__ failed.')

    # Sum
    if abs(e.get_sum()) > 0.0001:
        raise AssertionError('NiftiImageData3DDisplacement get_sum() failed.')

    # Dimensions
    f = e.get_dimensions()
    if not np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]):
        raise AssertionError('NiftiImageData3DDisplacement get_dimensions() failed.')

    # Get as array
    arr = d.as_array()
    if arr.max() != 100:
        raise AssertionError('NiftiImageData3DDisplacement as_array().max() failed.')
    if arr.ndim != 5:
        raise AssertionError('NiftiImageData3DDisplacement as_array() ndims failed.')
    if arr.shape != (64, 64, 64, 1, 3):
        raise AssertionError('NiftiImageData3DDisplacement as_array().shape failed.')

    # Check upsampling/downsampling
    u = pReg.NiftiImageData3DDisplacement(save_nifti_image_3d_displacement_not_split);
    original_spacing    = u.get_voxel_sizes();
    original_spacing    = original_spacing[1:4];
    upsampled_spacing   = [original_spacing[0]/2, original_spacing[1]/4, original_spacing[2]];
    downsampled_spacing = [original_spacing[0]*2, original_spacing[1]*4, original_spacing[2]];
    # Downsample
    v = u.deep_copy();
    v.set_voxel_spacing(downsampled_spacing,3);
    v.write(save_nifti_image_downsample);
    # Upsample then downsample, check nothing has changed
    w = u.deep_copy();
    w.set_voxel_spacing(upsampled_spacing,0);
    w.write(save_nifti_image_upsample);
    x = w.deep_copy();
    x.set_voxel_spacing(original_spacing,0);
    x.write(save_nifti_image_up_downsample);
    pReg.NiftiImageData.print_headers([u, v, w, x]);
    if x != u:
        raise AssertionError('NiftiImageData3DDisplacement::upsample()/downsample() failed.')

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImageData3DDisplacement test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImageData3DDeformation
def try_niftiimage3ddeformation():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData3DDeformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImageData3DDeformation from NiftiImageData3D
    b = pReg.NiftiImageData3DDeformation()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.write(save_nifti_image_3d_deformation_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_deformation_split)

    # Constructor from file
    c = pReg.NiftiImageData3DDeformation(save_nifti_image_3d_deformation_not_split)

    # Constructor from 3x3D
    d = pReg.NiftiImageData3DDeformation(ref_aladin, ref_aladin, ref_aladin)

    # Fill
    c.fill(100)

    # Get max
    if c.get_max() != 100:
        raise AssertionError('NiftiImageData3DDeformation fill()/get_max() failed.')

    # Get min
    if c.get_min() != 100:
        raise AssertionError('NiftiImageData3DDeformation fill()/get_min() failed.')

    # Deep copy
    d = c.deep_copy()
    if d.handle == c.handle:
        raise AssertionError('NiftiImageData3DDeformation deep_copy failed.')
    if d != c:
        raise AssertionError("NiftiImageData3DDeformation deep_copy failed.")

    # Addition
    e = d + d
    if abs(e.get_max() - 2 * d.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DDeformation __add__/get_max() failed.')

    # Subtraction
    e = d - d
    if abs(e.get_max()) > 0.0001:
        raise AssertionError('NiftiImageData3DDeformation __sub__ failed.')

    # Sum
    if abs(e.get_sum()) > 0.0001:
        raise AssertionError('NiftiImageData3DDeformation get_sum() failed.')

    # Dimensions
    f = e.get_dimensions()
    if not np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]):
        raise AssertionError('NiftiImageData3DDeformation get_dimensions() failed.')

    # Get as array
    arr = d.as_array()
    if arr.max() != 100:
        raise AssertionError('NiftiImageData3DDeformation as_array().max() failed.')
    if arr.ndim != 5:
        raise AssertionError('NiftiImageData3DDeformation as_array() ndims failed.')
    if arr.shape != (64, 64, 64, 1, 3):
        raise AssertionError('NiftiImageData3DDeformation as_array().shape failed.')

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImageData3DDeformation test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Nifty aladin
def try_niftyaladin():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty aladin test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # First set up some masks
    ref_mask = ref_aladin.deep_copy()
    flo_mask = flo_aladin.deep_copy()
    ref_mask.fill(1)
    flo_mask.fill(1)

    # default constructor
    na = pReg.NiftyAladinSym()
    na.set_reference_image(ref_aladin)
    na.set_floating_image(flo_aladin)
    na.set_parameter_file(parameter_file_aladin)
    na.set_parameter("SetInterpolationToCubic")
    na.set_parameter("SetLevelsToPerform", "1")
    na.set_parameter("SetMaxIterations", "5")
    na.set_reference_mask(ref_mask);
    na.set_floating_mask(flo_mask);
    na.process()

    # Get outputs
    warped = na.get_output()
    def_forward = na.get_deformation_field_forward()
    def_inverse = na.get_deformation_field_inverse()
    disp_forward = na.get_displacement_field_forward()
    disp_inverse = na.get_displacement_field_inverse()

    warped.write(aladin_warped)
    na.get_transformation_matrix_forward().write(TM_forward)
    na.get_transformation_matrix_inverse().write(TM_inverse)
    def_forward.write(aladin_def_forward)
    def_inverse.write_split_xyz_components(aladin_def_inverse)
    disp_forward.write(aladin_disp_forward)
    disp_inverse.write_split_xyz_components(aladin_disp_inverse)

    # forward TM
    forward_tm = na.get_transformation_matrix_forward()
    sys.stderr.write('\nforward tm:\n%s\n\n' % forward_tm.as_array())

    # Inverse TM
    inverse_tm = na.get_transformation_matrix_inverse()
    sys.stderr.write('\nInverse tm:\n%s\n\n' % inverse_tm.as_array())

    # Test converting disp to def
    a = pReg.NiftiImageData3DDeformation(disp_forward)
    if a != def_forward:
        raise AssertionError("NiftiImageData3DDeformation::create_from_disp() failed.")

    # Test converting def to disp
    b = pReg.NiftiImageData3DDisplacement(def_forward)
    if b != disp_forward:
        raise AssertionError("NiftiImageData3DDisplacement::create_from_def() failed.")

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty aladin test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    return na


# Nifty f3d
def try_niftyf3d():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty f3d test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Get initial transformation
    tm_init = pReg.AffineTransformation(TM_forward)

    # default constructor
    nf = pReg.NiftyF3dSym()
    nf.set_reference_image(ref_f3d)
    nf.set_floating_image(flo_f3d)
    nf.set_parameter_file(parameter_file_f3d)
    nf.set_reference_time_point(1)
    nf.set_floating_time_point(1)
    nf.set_initial_affine_transformation(tm_init)
    nf.process()

    # Get outputs
    warped = nf.get_output()
    def_forward = nf.get_deformation_field_forward()
    def_inverse = nf.get_deformation_field_inverse()
    disp_forward = nf.get_displacement_field_forward()
    disp_inverse = nf.get_displacement_field_inverse()

    warped.write(f3d_warped)
    def_forward.write(f3d_def_forward)
    def_inverse.write_split_xyz_components(f3d_def_inverse)
    disp_forward.write(f3d_disp_forward)
    disp_inverse.write_split_xyz_components(f3d_disp_inverse)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty f3d test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Transformation
def try_transformations(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Transformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Get transformations
    a3 = na.get_transformation_matrix_forward()
    b3 = na.get_displacement_field_forward()
    c3 = na.get_deformation_field_forward()

    # Get as deformations
    a_def = a3.get_as_deformation_field(ref_aladin)
    b_def = b3.get_as_deformation_field(ref_aladin)
    c_def = c3.get_as_deformation_field(ref_aladin)
    if a_def != na.get_deformation_field_forward():
        raise AssertionError()
    if b_def != na.get_deformation_field_forward():
        raise AssertionError()
    if c_def != na.get_deformation_field_forward():
        raise AssertionError()

    # Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = pReg.AffineTransformation.get_identity()
    trans = [tm_iden, tm_iden, c3]
    composed = pReg.NiftiImageData3DDeformation.compose_single_deformation(trans, ref_aladin)
    if composed != na.get_deformation_field_forward():
        raise AssertionError()

    # Test get_inverse
    tm_inv = tm_iden.get_inverse()

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Transformation test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Resample
def try_resample(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty resample test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    tm_iden = pReg.AffineTransformation.get_identity()
    tm      = na.get_transformation_matrix_forward()
    disp    = na.get_displacement_field_forward()
    deff    = na.get_deformation_field_forward()

    sys.stderr.write('Testing rigid resample...\n')
    nr1 = pReg.NiftyResample()
    nr1.set_reference_image(ref_aladin)
    nr1.set_floating_image(flo_aladin)
    nr1.set_interpolation_type_to_cubic_spline()  # try different interpolations
    nr1.set_interpolation_type(3)  # try different interpolations (cubic)
    nr1.add_transformation(tm_iden)
    nr1.add_transformation(tm)
    nr1.process()
    nr1.get_output().write(rigid_resample)

    sys.stderr.write('Testing non-rigid displacement...\n')
    nr2 = pReg.NiftyResample()
    nr2.set_reference_image(ref_aladin)
    nr2.set_floating_image(flo_aladin)
    nr2.set_interpolation_type_to_sinc()  # try different interpolations
    nr2.set_interpolation_type_to_linear()  # try different interpolations
    nr2.add_transformation(disp)
    nr2.process()
    nr2.get_output().write(nonrigid_resample_disp)

    sys.stderr.write('Testing non-rigid deformation...\n')
    nr3 = pReg.NiftyResample()
    nr3.set_reference_image(ref_aladin)
    nr3.set_floating_image(flo_aladin)
    nr3.set_interpolation_type_to_nearest_neighbour()  # try different interpolations
    nr3.add_transformation(deff)
    nr3.set_interpolation_type_to_linear()
    nr3.process()
    nr3.get_output().write(nonrigid_resample_def)

    # TODO this doesn't work. For some reason (even with NiftyReg directly), resampling with the TM from the registration
    # doesn't give the same result as the output from the registration itself (even with same interpolations). Even though 
    # ref and flo images are positive, the output of the registration can be negative. This implies that linear interpolation 
    # is not used to generate final image. You would hope it's used throughout the registration process, otherwise why is it there?
    # if na.get_output() != nr1.get_output():
    #     raise AssertionError()

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty resample test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Weighted mean
def try_weighted_mean(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting weighted mean test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Do 3D
    wm1 = pReg.ImageWeightedMean()
    # Change to float to avoid rounding errors
    im1 = ref_aladin.deep_copy()
    im2 = ref_aladin.deep_copy()
    im3 = ref_aladin.deep_copy()
    im4 = ref_aladin.deep_copy()
    im1.fill(1)
    im2.fill(4)
    im3.fill(7)
    im4.fill(6)
    wm1.add_image(im1, 2)
    wm1.add_image(im2, 4)
    wm1.add_image(im3, 3)
    wm1.add_image(im4, 1)
    wm1.process()
    wm1.get_output().write(output_weighted_mean)
    # Answer should be 4.5, so compare it to that!
    res = ref_aladin.deep_copy()
    res.fill(4.5)
    if wm1.get_output() != res:
        raise AssertionError()

    # Do 4D
    wm2 = pReg.ImageWeightedMean()
    im1 = na.get_deformation_field_forward().deep_copy()
    im2 = na.get_deformation_field_forward().deep_copy()
    im3 = na.get_deformation_field_forward().deep_copy()
    im4 = na.get_deformation_field_forward().deep_copy()
    im1.fill(1)
    im2.fill(4)
    im3.fill(7)
    im4.fill(6)
    wm2.add_image(im1, 2)
    wm2.add_image(im2, 4)
    wm2.add_image(im3, 3)
    wm2.add_image(im4, 1)
    wm2.process()
    wm2.get_output().write(output_weighted_mean_def)
    # Answer should be 4.5, so compare it to that!
    res = na.get_deformation_field_forward().deep_copy()
    res.fill(4.5)
    if wm2.get_output() != res:
        raise AssertionError()

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished weighted mean test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# AffineTransformation
def try_affinetransformation(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting AffineTransformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Construct from file
    a = pReg.AffineTransformation(TM_forward)
    if a.handle is None:
        raise AssertionError()

    # Multiply forward and inverse, should equal identity
    b = na.get_transformation_matrix_forward()
    c = na.get_transformation_matrix_inverse()
    d = b * c
    e = pReg.AffineTransformation.get_identity()
    if d != e:
        raise AssertionError('AffineTransformation::mult/comparison failed.')

    if e.get_determinant() - 1. > 1.e-7:
        raise AssertionError('AffineTransformation::get_determinant failed.')

    # Test get_Euler_angles
    array = np.zeros((4, 4), dtype=numpy.float32)

    array[0,2] =  1
    array[1,1] = -1
    array[2,0] = -1
    array[3,3] =  1
    test_Eul = pReg.AffineTransformation(array)
    # Example given by rotm2eul for MATLAB is [0 0 1; 0 -1 0; -1 0 0] -> XYZ = [-3.1416 1.5708 0]
    Eul = test_Eul.get_Euler_angles()
    Eul_expected = [-3.1416, 1.5708, 0]
    print(Eul)
    print(Eul_expected)
    if not np.allclose(Eul, Eul_expected, atol=1e-4):
        raise AssertionError('AffineTransformation get_Euler_angles() failed.')

    # Check as_array
    f = b.as_array()
    g = pReg.AffineTransformation(f)
    h = g.as_array()
    if not np.allclose(f, h, atol=1e-4):
        raise AssertionError('AffineTransformation as_array() failed.')

    # Average!
    trans = np.array([0., 0., 0.],dtype=numpy.float32)
    quat_1_array = np.array([0.92707, 0.02149, 0.19191, 0.32132],dtype=numpy.float32)
    quat_2_array = np.array([0.90361, 0.0025836, 0.097279, 0.41716],dtype=numpy.float32)
    quat_3_array = np.array([0.75868, -0.21289, 0.53263, 0.30884],dtype=numpy.float32)
    quat_1 = pReg.Quaternion(quat_1_array)
    quat_2 = pReg.Quaternion(quat_2_array)
    quat_3 = pReg.Quaternion(quat_3_array)
    tm_1 = pReg.AffineTransformation(trans,quat_1)
    tm_2 = pReg.AffineTransformation(trans,quat_2)
    tm_3 = pReg.AffineTransformation(trans,quat_3)
    average = pReg.AffineTransformation.get_average([tm_1, tm_2, tm_3])
    exptd_avg_array = np.zeros((4, 4), dtype=numpy.float32)
    exptd_avg_array[0][0] =  0.5836;
    exptd_avg_array[0][1] = -0.6736;
    exptd_avg_array[0][2] =  0.4535;
    exptd_avg_array[1][0] =  0.6007;
    exptd_avg_array[1][1] =  0.7339;
    exptd_avg_array[1][2] =  0.3171;
    exptd_avg_array[2][0] = -0.5464;
    exptd_avg_array[2][1] =  0.0874;
    exptd_avg_array[2][2] =  0.8329;
    exptd_avg_array[3][3] =  1;
    exptd_average = pReg.AffineTransformation(exptd_avg_array)
    if exptd_average != average, atol=1e-4:
        raise AssertionError('AffineTransformation average failed.')
    print(average.as_array())


    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished AffineTransformation test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

# Quaternion
def try_quaternion():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Quaternion test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Construct TM
    array = np.zeros((4, 4), dtype=numpy.float32)
    array[0,2] = 1
    array[1,1] = 1
    array[2,0] = -1
    array[3,3] = 1
    rotm = pReg.AffineTransformation(array)

    # Convert to quaternion
    quat = pReg.Quaternion(rotm)
    a = quat.as_array()

    # Construct from numpy array
    expt_array = np.array([0.707107, 0., 0.707107, 0.],dtype=numpy.float32)
    expt = pReg.Quaternion(expt_array)

    # Compare to expected values
    if not np.allclose(quat.as_array(), expt_array, atol=1e-4):
        raise AssertionError('Quaternion from TM failed.')

    # Convert back to TM
    trans_array = np.array([0., 0., 0.],dtype=numpy.float32)
    affine = pReg.AffineTransformation(trans_array,quat)
    if affine != rotm:
        raise AssertionError('TM to quaternion failed.')

    # Convert TM to quaternion
    quat2 = affine.get_quaternion()
    if not np.allclose(quat.as_array(), quat2.as_array(), atol=1e-4):
        raise AssertionError('AffineTransformation:get_quaternion() failed.')

    # Average!
    quat_1_array = np.array([0.92707, 0.02149, 0.19191, 0.32132],dtype=numpy.float32)
    quat_2_array = np.array([0.90361, 0.0025836, 0.097279, 0.41716],dtype=numpy.float32)
    quat_3_array = np.array([0.75868, -0.21289, 0.53263, 0.30884],dtype=numpy.float32)
    quat_1 = pReg.Quaternion(quat_1_array)
    quat_2 = pReg.Quaternion(quat_2_array)
    quat_3 = pReg.Quaternion(quat_3_array)
    exptd_avg_array = np.array([0.88748, -0.0647152, 0.281671, 0.35896],dtype=numpy.float32)
    exptd_average = pReg.Quaternion(exptd_avg_array)
    average = pReg.Quaternion.get_average([quat_1, quat_2, quat_3])
    if not np.allclose(exptd_average.as_array(), average.as_array(), atol=1e-4):
        raise AssertionError('Quaternion average failed.')
    print(average.as_array())


    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Quaternion test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
    try_niftiimage()
    try_niftiimage3d()
    try_niftiimage3dtensor()
    try_niftiimage3ddisplacement()
    try_niftiimage3ddeformation()
    na = try_niftyaladin()
    try_niftyf3d()
    try_transformations(na)
    try_resample(na)
    try_weighted_mean(na)
    try_affinetransformation(na)
    try_quaternion()


if __name__ == "__main__":
    try:
        test()
    except:
        raise error("Error encountered.")
