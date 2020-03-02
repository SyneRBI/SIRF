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
import nibabel as nib
import sirf.Reg
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
aladin_def_inverse_xyz = output_prefix + "aladin_def_inverse_%s.nii"
aladin_def_inverse = output_prefix + "aladin_def_inverse.nii"
aladin_def_fwd_inv = output_prefix + "aladin_def_fwd_then_inv.nii"
aladin_disp_forward = output_prefix + "aladin_disp_forward.nii"
aladin_disp_inverse = output_prefix + "aladin_disp_inverse_%s.nii"
f3d_def_forward = output_prefix + "f3d_disp_forward.nii"
f3d_def_inverse = output_prefix + "f3d_disp_inverse_%s.nii"
f3d_disp_forward = output_prefix + "f3d_disp_forward.nii"
f3d_disp_inverse = output_prefix + "f3d_disp_inverse_%s.nii"

rigid_resample = output_prefix + "rigid_resample.nii"
nonrigid_resample_disp = output_prefix + "nonrigid_resample_disp.nii"
nonrigid_resample_def = output_prefix + "nonrigid_resample_def.nii"
niftymomo_resample_adj = output_prefix + "niftymomo_resample_adj.nii"
output_weighted_mean = output_prefix + "weighted_mean.nii"
output_weighted_mean_def = output_prefix + "weighted_mean_def.nii"
output_float = output_prefix + "reg_aladin_float.nii"

ref_aladin = sirf.Reg.NiftiImageData3D(ref_aladin_filename)
flo_aladin = sirf.Reg.NiftiImageData3D(flo_aladin_filename)
ref_f3d = sirf.Reg.NiftiImageData3D(ref_f3d_filename)
flo_f3d = sirf.Reg.NiftiImageData3D(flo_f3d_filename)

# NiftiImageData
def try_niftiimage():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = sirf.Reg.NiftiImageData()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = sirf.Reg.NiftiImageData(ref_aladin_filename)

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
    ref_aladin_float = sirf.Reg.NiftiImageData3D(output_float)
    arr1 = ref_aladin.as_array()
    arr2 = ref_aladin_float.as_array()
    if not np.array_equal(arr1,arr2):
        raise AssertionError("NiftiImageData::write()/change_datatype() failed.")

    # Test print methods
    q.print_header()
    sirf.Reg.NiftiImageData.print_headers([q, s])

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
    u = sirf.Reg.NiftiImageData(ref_aladin_filename);
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
    sirf.Reg.NiftiImageData.print_headers([u, v, w, x]);
    if x != u:
        raise AssertionError('NiftiImageData::upsample()/downsample() failed.')

    # Check get_contains_nans
    x_arr = x.as_array()
    x_arr.fill(0)
    x.fill(x_arr)
    if x.get_contains_nans():
        raise AssertionError('NiftiImageData::get_contains_nans() 1 failed.')
    x_arr[1] = np.nan
    x.fill(x_arr)
    if not x.get_contains_nans():
        raise AssertionError('NiftiImageData::get_contains_nans() 2 failed.')

    # Test that fill works regardless of C- or F-style numpy arrays
    im = sirf.Reg.NiftiImageData(ref_aladin_filename)
    arr = im.as_array()
    arr_C = numpy.ascontiguousarray(arr)
    arr_F = numpy.asfortranarray(arr)
    im.fill(arr_C)
    arr_C2 = im.as_array()
    im.fill(arr_F)
    arr_F2 = im.as_array()
    if not np.array_equal(arr_C2, arr_F2):
        raise AssertionError('NiftiImageData::fill() failed for C- or F-style numpy arrays.')

    # Compare between sirf.Reg.NiftiImageData::as_array() and nibabel
    arr1 = sirf.Reg.NiftiImageData(ref_aladin_filename).as_array()
    arr2 = nib.load(ref_aladin_filename).get_fdata()
    if not numpy.array_equal(arr1,arr2):
        raise AssertionError("NiftiImageData as_array() failed.")

    # Test geom info
    geom_info = im.get_geometrical_info()
    geom_info.print_info()
    # Get voxel sizes
    if geom_info.get_size() != (64, 64, 64):
        raise AssertionError("SIRF get_geometrical_info().get_size() failed.")
    if geom_info.get_spacing() != (4.0625, 4.0625, 4.0625):
        raise AssertionError("SIRF get_geometrical_info().get_spacing() failed.")

    im.standardise()
    if abs(im.get_standard_deviation() - 1) > 0.01:
        raise AssertionError("NiftiImageData standardise() or get_standard_deviation() failed.")
    if abs(im.get_variance() - 1) > 0.01:
        raise AssertionError("NiftiImageData standardise() or get_variance() failed.")
    if abs(im.get_mean()) > 0.0001:
        raise AssertionError("NiftiImageData standardise() or get_mean() failed.")

    # Check normalise 
    im.normalise_zero_and_one()
    if abs(im.get_min()) > 0.0001 or abs(im.get_max()-1) > 0.0001:
        raise AssertionError("NiftiImageData normalise_between_zero_and_one() failed.")

    # Test inner product
    in1 = x.deep_copy()
    in2 = x.deep_copy()
    in1_arr = in1.as_array()
    in2_arr = in2.as_array()
    dims = in1.get_dimensions()
    inner_product = 0
    for idx_x in range(dims[1]):
        for idx_y in range(dims[2]):
            for idx_z in range(dims[3]):
                in1_arr[idx_x, idx_y, idx_z] = float(i)
                in2_arr[idx_x, idx_y, idx_z] = float(3*i-1)
                inner_product += float(i) * float(3*i-1)
    in1.fill(in1_arr)
    in2.fill(in2_arr)
    if abs(inner_product - in1.get_inner_product(in2)) > 1e-4:
        raise AssertionError("NiftiImageData::get_inner_product() failed.")

    # Pad then crop, should be the same
    aa = ref_aladin
    cc = aa.clone()
    original_dims = aa.get_dimensions()

    pad_in_min_dir = [1, 2, 3, 0, 0, 0, 0]
    pad_in_max_dir = [4, 5, 6, 0, 0, 0, 0]
    cc.pad(pad_in_min_dir, pad_in_max_dir, 100.)

    padded_dims = cc.get_dimensions()
    for i in range(7):
        if padded_dims[i+1] != original_dims[i+1] + pad_in_min_dir[i] + pad_in_max_dir[i]:
            raise AssertionError("NiftiImageData::pad failed")

    # Crop back to beginning
    cropped_min_dir = pad_in_min_dir
    cropped_max_dir = list(cropped_min_dir)
    for i in range(7):
        cropped_max_dir[i] = original_dims[i+1] + cropped_min_dir[i] - 1

    cc.crop(cropped_min_dir, cropped_max_dir)
    if aa != cc:
        raise AssertionError("NiftiImageData::pad/crop failed")

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
    a = sirf.Reg.NiftiImageData3D()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = sirf.Reg.NiftiImageData3D(ref_aladin_filename)

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

    # try linear algebra
    h = d/10000;
    if abs(h.get_max()-d.get_max()/10000) > 1e-4:
        raise AssertionError('NiftiImageData3D linear algebra failed.')

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
    b = sirf.Reg.NiftiImageData3DTensor()
    b.create_from_3D_image(ref_aladin)

    # # Save to file
    b.write(save_nifti_image_3d_tensor_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_tensor_split)

    # Constructor from file
    c = sirf.Reg.NiftiImageData3DTensor(save_nifti_image_3d_tensor_not_split)

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
    h = sirf.Reg.NiftiImageData3DTensor(im1, im2, im3)

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
    b = sirf.Reg.NiftiImageData3DDisplacement()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.write(save_nifti_image_3d_displacement_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_displacement_split)

    # Constructor from file
    c = sirf.Reg.NiftiImageData3DDisplacement(save_nifti_image_3d_displacement_not_split)

    # Constructor from 3x3D
    d = sirf.Reg.NiftiImageData3DDisplacement(ref_aladin, ref_aladin, ref_aladin)

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
    u = sirf.Reg.NiftiImageData3DDisplacement(save_nifti_image_3d_displacement_not_split);
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
    sirf.Reg.NiftiImageData.print_headers([u, v, w, x]);
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
    b = sirf.Reg.NiftiImageData3DDeformation()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.write(save_nifti_image_3d_deformation_not_split)
    b.write_split_xyz_components(save_nifti_image_3d_deformation_split)

    # Constructor from file
    c = sirf.Reg.NiftiImageData3DDeformation(save_nifti_image_3d_deformation_not_split)

    # Constructor from 3x3D
    d = sirf.Reg.NiftiImageData3DDeformation(ref_aladin, ref_aladin, ref_aladin)

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

    # Print all wrapped methods.
    sirf.Reg.NiftyAladinSym.print_all_wrapped_methods()

    # default constructor
    na = sirf.Reg.NiftyAladinSym()
    na.set_reference_image(ref_aladin)
    na.set_floating_image(flo_aladin)
    na.set_parameter_file(parameter_file_aladin)
    na.set_parameter("SetInterpolationToCubic")
    na.set_parameter("SetLevelsToPerform", "1")
    na.set_parameter("SetMaxIterations", "5")
    na.set_parameter("SetPerformRigid", "1")
    na.set_parameter("SetPerformAffine", "0")
    na.set_reference_mask(ref_mask)
    na.set_floating_mask(flo_mask)
    na.process()

    # Get outputs
    warped = na.get_output().deep_copy()
    def_forward = na.get_deformation_field_forward().deep_copy()
    def_inverse = na.get_deformation_field_inverse().deep_copy()
    disp_forward = na.get_displacement_field_forward().deep_copy()
    disp_inverse = na.get_displacement_field_inverse().deep_copy()
    TM_forward_ = na.get_transformation_matrix_forward().deep_copy()
    TM_inverse_ = na.get_transformation_matrix_inverse().deep_copy()

    # Test via filenames
    na.set_reference_image_filename(ref_aladin_filename)
    na.set_floating_image_filename(flo_aladin_filename)
    na.process()

    if warped != na.get_output() or \
        def_forward != na.get_deformation_field_forward() or \
        def_inverse != na.get_deformation_field_inverse() or \
        disp_forward != na.get_displacement_field_forward() or \
        disp_inverse != na.get_displacement_field_inverse() or \
        TM_forward_ != na.get_transformation_matrix_forward() or \
        TM_inverse_ != na.get_transformation_matrix_inverse():
        raise AssertionError()

    warped.write(aladin_warped)
    TM_forward_.write(TM_forward)
    TM_inverse_.write(TM_inverse)
    def_forward.write(aladin_def_forward)
    def_inverse.write_split_xyz_components(aladin_def_inverse_xyz)
    def_inverse.write(aladin_def_inverse)
    disp_forward.write(aladin_disp_forward)
    disp_inverse.write_split_xyz_components(aladin_disp_inverse)

    # forward TM
    forward_tm = na.get_transformation_matrix_forward()
    sys.stderr.write('\nforward tm:\n%s\n\n' % forward_tm.as_array())

    # Inverse TM
    inverse_tm = na.get_transformation_matrix_inverse()
    sys.stderr.write('\nInverse tm:\n%s\n\n' % inverse_tm.as_array())

    # Test converting disp to def
    a = sirf.Reg.NiftiImageData3DDeformation(disp_forward)
    if a != def_forward:
        raise AssertionError("NiftiImageData3DDeformation::create_from_disp() failed.")

    # Test converting def to disp
    b = sirf.Reg.NiftiImageData3DDisplacement(def_forward)
    if b != disp_forward:
        raise AssertionError("NiftiImageData3DDisplacement::create_from_def() failed.")

    # Check NiftiImageData3DDeformation::get_inverse()
    def_fwd_then_inv = def_forward.get_inverse(flo_aladin)
    def_fwd_then_inv.write(aladin_def_fwd_inv)
    sirf.Reg.NiftiImageData.print_headers([ref_aladin, flo_aladin, def_inverse, def_fwd_then_inv])

    # Reference forward with def_inv
    resample = sirf.Reg.NiftyResample()
    resample.set_reference_image(flo_aladin)
    resample.set_floating_image(ref_aladin)
    resample.set_padding_value(0.)
    resample.set_interpolation_type_to_linear()
    resample.add_transformation(def_inverse)
    out1 = resample.forward(ref_aladin)

    # Reference forward with def_fwd_then_inv
    resample.clear_transformations()
    resample.add_transformation(def_fwd_then_inv)
    out2 = resample.forward(ref_aladin)

    sirf.Reg.NiftiImageData.print_headers([out1, out2])
    if out1 != out2:
        raise AssertionError("NiftiImageData3DDeformation::get_inverse() failed.")

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
    tm_init = sirf.Reg.AffineTransformation(TM_forward)

    # Print all wrapped methods.
    sirf.Reg.NiftyF3dSym.print_all_wrapped_methods()

    # default constructor
    nf = sirf.Reg.NiftyF3dSym()
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
    tm_iden = sirf.Reg.AffineTransformation.get_identity()
    trans = [tm_iden, tm_iden, c3]
    composed = sirf.Reg.NiftiImageData3DDeformation.compose_single_deformation(trans, ref_aladin)
    if composed != na.get_deformation_field_forward():
        raise AssertionError()

    # Test get_inverse
    tm_inv = tm_iden.get_inverse()
    if not tm_inv:
        raise AssertionError()

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

    tm_iden = sirf.Reg.AffineTransformation.get_identity()
    tm      = na.get_transformation_matrix_forward()
    disp    = na.get_displacement_field_forward()
    deff    = na.get_deformation_field_forward()
    padding_value = -20

    sys.stderr.write('Testing rigid resample...\n')
    nr1 = sirf.Reg.NiftyResample()
    nr1.set_reference_image(ref_aladin)
    nr1.set_floating_image(flo_aladin)
    nr1.set_interpolation_type_to_cubic_spline()  # try different interpolations
    nr1.set_interpolation_type(3)  # try different interpolations (cubic)
    nr1.add_transformation(tm_iden)
    nr1.clear_transformations()
    nr1.add_transformation(tm_iden)
    nr1.add_transformation(tm)
    nr1.process()
    nr1.get_output().write(rigid_resample)

    sys.stderr.write('Testing non-rigid displacement...\n')
    nr2 = sirf.Reg.NiftyResample()
    nr2.set_reference_image(ref_aladin)
    nr2.set_floating_image(flo_aladin)
    nr2.set_interpolation_type_to_sinc()  # try different interpolations
    nr2.set_interpolation_type_to_nearest_neighbour()  # try different interpolations
    nr2.add_transformation(disp)
    nr2.set_padding_value(padding_value)
    nr2.process()
    nr2.get_output().write(nonrigid_resample_disp)

    if nr2.get_output().get_min() != padding_value:
        raise AssertionError('NiftyResample:set_padding_value failed')

    sys.stderr.write('Testing non-rigid deformation...\n')
    nr3 = sirf.Reg.NiftyResample()
    nr3.set_reference_image(ref_aladin)
    nr3.set_floating_image(flo_aladin)
    nr3.set_interpolation_type_to_linear()  # try different interpolations
    nr3.add_transformation(deff)
    nr3.set_interpolation_type_to_linear()
    nr3.process()
    nr3.get_output().write(nonrigid_resample_def)

    # Check that the following give the same result
    #       out = resample.forward(in)
    #       resample.forward(out, in)
    out1 = nr3.forward(flo_aladin)
    out2 = ref_aladin.deep_copy()
    nr3.forward(out2, flo_aladin)
    if out1 != out2:
        raise AssertionError('out = NiftyResample::forward(in) and NiftyResample::forward(out, in) do not give same result.')

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


# NiftyMoMo
def try_niftymomo(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftyMomMo test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # The forward and the adjoint should meet the following criterion:
    # | < x, Ty > - < y, Tsx > | / 0.5 * (| < x, Ty > | + | < y, Tsx > |) < epsilon
    # for all images x and y, where T is the transform and Ts is the adjoint.

    x = ref_aladin
    T = na.get_transformation_matrix_forward()
    y = flo_aladin

    # Add in a magnification to make things interesting
    t = T.as_array()
    t[0][0] = 1.5
    T = sirf.Reg.AffineTransformation(t)

    # make it slightly unsquare to spice things up
    min_idx = [0, 1, 2]
    y_dims = y.get_dimensions()
    max_idx = [y_dims[1] - 3, y_dims[2] - 1, y_dims[3]-5]
    y.crop(min_idx, max_idx)

    sys.stderr.write('Testing adjoint resample...\n')
    nr = sirf.Reg.NiftyResample()
    nr.set_reference_image(x)
    nr.set_floating_image(y)
    nr.set_interpolation_type_to_linear()
    nr.add_transformation(T)

    # Do the forward
    Ty = nr.forward(y)

    # Do the adjoint
    Tsx = nr.adjoint(x)

    # Check the adjoint is truly the adjoint with: |<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) < epsilon
    inner_x_Ty = x.get_inner_product(Ty)
    inner_y_Tsx = y.get_inner_product(Tsx)
    adjoint_test = abs(inner_x_Ty - inner_y_Tsx) / (0.5 * (abs(inner_x_Ty) + abs(inner_y_Tsx)))
    sys.stderr.write('<x, Ty>  = %f\n' % inner_x_Ty)
    sys.stderr.write('<y, Tsx> = %f\n' % inner_y_Tsx)
    sys.stderr.write('|<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) = %f\n' % adjoint_test)
    if adjoint_test > 1e-4:
        raise AssertionError("NiftyResample::adjoint() failed")

    # Check that the following give the same result
    #       out = resample.adjoint(in)
    #       resample.adjoint(out, in)
    out1 = nr.adjoint(x)
    out2 = y.deep_copy()
    nr.backward(out2, x)
    if out1 != out2:
        raise AssertionError(
            'out = NiftyResample::adjoint(in) and NiftyResample::adjoint(out, in) do not give same result.')

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftyMoMo test.\n')
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
    wm1 = sirf.Reg.ImageWeightedMean()
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
    wm2 = sirf.Reg.ImageWeightedMean()
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
    a = sirf.Reg.AffineTransformation(TM_forward)
    if a.handle is None:
        raise AssertionError()

    # Multiply forward and inverse, should equal identity
    b = na.get_transformation_matrix_forward()
    c = na.get_transformation_matrix_inverse()
    d = b * c
    e = sirf.Reg.AffineTransformation.get_identity()
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
    test_Eul = sirf.Reg.AffineTransformation(array)
    # Example given by rotm2eul for MATLAB is [0 0 1; 0 -1 0; -1 0 0] -> XYZ = [-3.1416 1.5708 0]
    Eul = test_Eul.get_Euler_angles()
    Eul_expected = [-3.1416, 1.5708, 0]
    print(Eul)
    print(Eul_expected)
    if not np.allclose(Eul, Eul_expected, atol=1e-4):
        raise AssertionError('AffineTransformation get_Euler_angles() failed.')

    # Check as_array
    f = b.as_array()
    g = sirf.Reg.AffineTransformation(f)
    h = g.as_array()
    if not np.allclose(f, h, atol=1e-4):
        raise AssertionError('AffineTransformation as_array() failed.')

    # Average!
    trans = np.array([0., 0., 0.],dtype=numpy.float32)
    quat_1_array = np.array([0.92707, 0.02149, 0.19191, 0.32132],dtype=numpy.float32)
    quat_2_array = np.array([0.90361, 0.0025836, 0.097279, 0.41716],dtype=numpy.float32)
    quat_3_array = np.array([0.75868, -0.21289, 0.53263, 0.30884],dtype=numpy.float32)
    quat_1 = sirf.Reg.Quaternion(quat_1_array)
    quat_2 = sirf.Reg.Quaternion(quat_2_array)
    quat_3 = sirf.Reg.Quaternion(quat_3_array)
    tm_1 = sirf.Reg.AffineTransformation(trans,quat_1)
    tm_2 = sirf.Reg.AffineTransformation(trans,quat_2)
    tm_3 = sirf.Reg.AffineTransformation(trans,quat_3)
    average = sirf.Reg.AffineTransformation.get_average([tm_1, tm_2, tm_3])
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
    exptd_average = sirf.Reg.AffineTransformation(exptd_avg_array)
    if exptd_average != average:
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
    rotm = sirf.Reg.AffineTransformation(array)

    # Convert to quaternion
    quat = sirf.Reg.Quaternion(rotm)
    a = quat.as_array()
    if a is None:
        raise AssertionError()

    # Construct from numpy array
    expt_array = np.array([0.707107, 0., 0.707107, 0.],dtype=numpy.float32)
    expt = sirf.Reg.Quaternion(expt_array)
    if expt is None:
        raise AssertionError()

    # Compare to expected values
    if not np.allclose(quat.as_array(), expt_array, atol=1e-4):
        raise AssertionError('Quaternion from TM failed.')

    # Convert back to TM
    trans_array = np.array([0., 0., 0.],dtype=numpy.float32)
    affine = sirf.Reg.AffineTransformation(trans_array,quat)
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
    quat_1 = sirf.Reg.Quaternion(quat_1_array)
    quat_2 = sirf.Reg.Quaternion(quat_2_array)
    quat_3 = sirf.Reg.Quaternion(quat_3_array)
    exptd_avg_array = np.array([0.88748, -0.0647152, 0.281671, 0.35896],dtype=numpy.float32)
    exptd_average = sirf.Reg.Quaternion(exptd_avg_array)
    average = sirf.Reg.Quaternion.get_average([quat_1, quat_2, quat_3])
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
    try_niftymomo(na)
    try_weighted_mean(na)
    try_affinetransformation(na)
    try_quaternion()


if __name__ == "__main__":
    try:
        test()
    except:
        raise error("Error encountered.")
