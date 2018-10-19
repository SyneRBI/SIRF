# Imports
import os
import sys
import time

import numpy as np
import pSIRFReg
import pSTIR
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

ref_aladin = pSIRFReg.NiftiImageData3D(ref_aladin_filename)
flo_aladin = pSIRFReg.NiftiImageData3D(flo_aladin_filename)
ref_f3d = pSIRFReg.NiftiImageData3D(ref_f3d_filename)
flo_f3d = pSIRFReg.NiftiImageData3D(flo_f3d_filename)

# NiftiImageData
def try_niftiimage():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImageData test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = pSIRFReg.NiftiImageData()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = pSIRFReg.NiftiImageData(ref_aladin_filename)

    # Save to file
    b.save_to_file(save_nifti_image)

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
    ref_aladin.save_to_file(output_float, 16) # float
    ref_aladin_float = pSIRFReg.NiftiImageData3D(output_float)
    arr1 = ref_aladin.as_array()
    arr2 = ref_aladin_float.as_array()
    if not np.array_equal(arr1,arr2):
        raise AssertionError("SIRFRegMisc::save_to_file()/change_datatype() failed.")

    # Test print methods
    q.print_header()
    pSIRFReg.NiftiImageData.print_headers([q, s])

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
    a = pSIRFReg.NiftiImageData3D()
    if a.handle is None:
        raise AssertionError()

    # Read from file
    b = pSIRFReg.NiftiImageData3D(ref_aladin_filename)

    # Save to file
    b.save_to_file(save_nifti_image_3d)

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
    b = pSIRFReg.NiftiImageData3DTensor()
    b.create_from_3D_image(ref_aladin)

    # # Save to file
    b.save_to_file(save_nifti_image_3d_tensor_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_tensor_split)

    # Constructor from file
    c = pSIRFReg.NiftiImageData3DTensor(save_nifti_image_3d_tensor_not_split)

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
    h = pSIRFReg.NiftiImageData3DTensor(im1, im2, im3)

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
    b = pSIRFReg.NiftiImageData3DDisplacement()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.save_to_file(save_nifti_image_3d_displacement_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_displacement_split)

    # Constructor from file
    c = pSIRFReg.NiftiImageData3DDisplacement(save_nifti_image_3d_displacement_not_split)

    # Constructor from 3x3D
    d = pSIRFReg.NiftiImageData3DDisplacement(ref_aladin, ref_aladin, ref_aladin)

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
    b = pSIRFReg.NiftiImageData3DDeformation()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.save_to_file(save_nifti_image_3d_deformation_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_deformation_split)

    # Constructor from file
    c = pSIRFReg.NiftiImageData3DDeformation(save_nifti_image_3d_deformation_not_split)

    # Constructor from 3x3D
    d = pSIRFReg.NiftiImageData3DDeformation(ref_aladin, ref_aladin, ref_aladin)

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

    # default constructor
    na = pSIRFReg.NiftyAladinSym()
    na.set_reference_image(ref_aladin)
    na.set_floating_image(flo_aladin)
    na.set_parameter_file(parameter_file_aladin)
    na.set_parameter("SetInterpolationToCubic")
    na.set_parameter("SetLevelsToPerform", "1")
    na.set_parameter("SetMaxIterations", "5")
    na.process()

    # Get outputs
    warped = na.get_output()
    def_forward = na.get_deformation_field_forward()
    def_inverse = na.get_deformation_field_inverse()
    disp_forward = na.get_displacement_field_forward()
    disp_inverse = na.get_displacement_field_inverse()

    warped.save_to_file(aladin_warped)
    na.get_transformation_matrix_forward().save_to_file(TM_forward)
    na.get_transformation_matrix_inverse().save_to_file(TM_inverse)
    def_forward.save_to_file(aladin_def_forward)
    def_inverse.save_to_file_split_xyz_components(aladin_def_inverse)
    disp_forward.save_to_file(aladin_disp_forward)
    disp_inverse.save_to_file_split_xyz_components(aladin_disp_inverse)

    # forward TM
    forward_tm = na.get_transformation_matrix_forward()
    sys.stderr.write('\nforward tm:\n%s\n\n' % forward_tm.as_array())

    # Inverse TM
    inverse_tm = na.get_transformation_matrix_inverse()
    sys.stderr.write('\nInverse tm:\n%s\n\n' % inverse_tm.as_array())

    # Test converting disp to def
    a = pSIRFReg.NiftiImageData3DDeformation()
    a.create_from_disp(disp_forward)
    if a != def_forward:
        raise AssertionError("NiftiImageData3DDeformation::create_from_disp() failed.")

    # Test converting def to disp
    b = pSIRFReg.NiftiImageData3DDisplacement()
    b.create_from_def(def_forward)
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
    tm_init = pSIRFReg.AffineTransformation(TM_forward)

    # default constructor
    nf = pSIRFReg.NiftyF3dSym()
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

    warped.save_to_file(f3d_warped)
    def_forward.save_to_file(f3d_def_forward)
    def_inverse.save_to_file_split_xyz_components(f3d_def_inverse)
    disp_forward.save_to_file(f3d_disp_forward)
    disp_inverse.save_to_file_split_xyz_components(f3d_disp_inverse)

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
    tm_iden = pSIRFReg.AffineTransformation.get_identity()
    trans = [tm_iden, tm_iden, c3]
    composed = pSIRFReg.NiftiImageData3DDeformation.compose_single_deformation(trans, ref_aladin)
    if composed != na.get_deformation_field_forward():
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

    tm_iden = pSIRFReg.AffineTransformation.get_identity()
    tm      = na.get_transformation_matrix_forward()
    disp    = na.get_displacement_field_forward()
    deff    = na.get_deformation_field_forward()

    sys.stderr.write('Testing rigid resample...\n')
    nr1 = pSIRFReg.NiftyResample()
    nr1.set_reference_image(ref_aladin)
    nr1.set_floating_image(flo_aladin)
    nr1.set_interpolation_type_to_cubic_spline()  # try different interpolations
    nr1.set_interpolation_type(3)  # try different interpolations (cubic)
    nr1.add_transformation_affine(tm_iden)
    nr1.add_transformation_affine(tm)
    nr1.process()
    nr1.get_output().save_to_file(rigid_resample)

    sys.stderr.write('Testing non-rigid displacement...\n')
    nr2 = pSIRFReg.NiftyResample()
    nr2.set_reference_image(ref_aladin)
    nr2.set_floating_image(flo_aladin)
    nr2.set_interpolation_type_to_sinc()  # try different interpolations
    nr2.set_interpolation_type_to_linear()  # try different interpolations
    nr2.add_transformation_disp(disp)
    nr2.process()
    nr2.get_output().save_to_file(nonrigid_resample_disp)

    sys.stderr.write('Testing non-rigid deformation...\n')
    nr3 = pSIRFReg.NiftyResample()
    nr3.set_reference_image(ref_aladin)
    nr3.set_floating_image(flo_aladin)
    nr3.set_interpolation_type_to_nearest_neighbour()  # try different interpolations
    nr3.add_transformation_def(deff)
    nr3.set_interpolation_type_to_linear()
    nr3.process()
    nr3.get_output().save_to_file(nonrigid_resample_def)

    if na.get_output() != nr1.get_output():
        raise AssertionError()

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
    wm1 = pSIRFReg.ImageWeightedMean()
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
    wm1.get_output().save_to_file(output_weighted_mean)
    # Answer should be 4.5, so compare it to that!
    res = ref_aladin.deep_copy()
    res.fill(4.5)
    if wm1.get_output() != res:
        raise AssertionError()

    # Do 4D
    wm2 = pSIRFReg.ImageWeightedMean()
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
    wm2.get_output().save_to_file(output_weighted_mean_def)
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


# STIR to SIRFReg
def try_stir_to_sirfreg():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting STIR to SIRFReg test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Open stir image
    pet_image_data = pSTIR.ImageData(ref_aladin_filename)
    image_data_from_stir = pSIRFReg.NiftiImageData3D(pet_image_data)

    # Now fill the stir and sirfreg images with 1 and 100, respectively
    pet_image_data.fill(1.)
    image_data_from_stir.fill(100)
    if pet_image_data.as_array().max() == image_data_from_stir.get_max():
        raise AssertionError()

    # Fill the stir image with the sirfreg
    image_data_from_stir.copy_data_to(pet_image_data)
    if pet_image_data.as_array().max() != image_data_from_stir.get_max():
        raise AssertionError()

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished STIR to SIRFReg test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# SIRFRegAffineTransformation
def try_sirfregAffineTransformation(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting SIRFRegAffineTransformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Construct from file
    a = pSIRFReg.AffineTransformation(TM_forward)
    if a.handle is None:
        raise AssertionError()

    # Multiply forward and inverse, should equal identity
    b = na.get_transformation_matrix_forward()
    c = na.get_transformation_matrix_inverse()
    d = b * c
    e = pSIRFReg.AffineTransformation.get_identity()
    if d != e:
        raise AssertionError('SIRFRegAffineTransformation::mult/comparison failed.')

    if e.get_determinant() - 1. > 1.e-7:
        raise AssertionError('SIRFRegAffineTransformation::get_determinant failed.')

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished SIRFRegAffineTransformation test.\n')
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
    try_stir_to_sirfreg()
    try_sirfregAffineTransformation(na)


if __name__ == "__main__":
    test()
