# -*- coding: utf-8 -*-
"""sirf.Gadgetron Test set 1.
v{version}

Fully sampled data tests

Usage:
  test1 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import numpy as np
from numpy.core.function_base import linspace
from pathlib import Path

import sirf.DynamicSimulation as pDS
import sirf.Gadgetron as pMR
import sirf.Reg as pReg

from sirf.Utilities import __license__, runner

__version__ = "3.1.0"
__author__ = "Johannes Mayer"



def prepare_test_simulation(fname_mr_rawdata, fpath_xml):

    rawdata = pMR.AcquisitionData(fname_mr_rawdata)
    rawdata = pMR.preprocess_acquisition_data(rawdata)

    empty_img = pMR.ImageData()
    empty_img.from_acquisition_data (rawdata)

    labels = pReg.NiftiImageData3D(empty_img)

    label_array = labels.as_array()
    Nx, Ny, Nz = label_array.shape

    tissue_label = 13

    label_array[:] = 0
    label_array[Nx//4:3*Nx//4, Ny//4:3*Ny//4, Nz//4:3*Nz//4]  = tissue_label

    labels.fill(label_array)

    mrsim = pDS.DynamicSimulation(labels, fpath_xml)
    
    return mrsim, rawdata, labels


def prep_displacement_field(nifti_3D_volume):
    
    dvf_x = nifti_3D_volume.deep_copy()
    dvf_y = nifti_3D_volume.deep_copy()
    dvf_z = nifti_3D_volume.deep_copy()

    Nx, Ny, Nz = nifti_3D_volume.shape

    dvf_x.fill(Nx / 3)
    dvf_y.fill(Ny / 3)
    dvf_z.fill(Nz / 3)

    return pReg.NiftiImageData3DDisplacement(dvf_x, dvf_y, dvf_z)


def test_static_mr_simulation(rec=False, verb=False, throw=True):

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/cDynamicSimulation/'

    fpath_xml = input_fpath_prefix + 'Segmentations/XCAT_TissueParameters_XML.xml'
    fpath_template_rawdata = input_fpath_prefix + 'TemplateData/MR/CV_nav_cart_64Cube_1Echo.h5'

    mrsim, rawdata, __ = prepare_test_simulation(fpath_template_rawdata, fpath_xml)
    
    mrsim.set_acquisition_template_data(rawdata)

    csm = pMR.CoilSensitivityData()
    csm.calculate(rawdata)
    mrsim.set_csm(csm)

    SNR = 5
    SNR_label = 13

    mrsim.set_snr(SNR)
    mrsim.set_snr_label(SNR_label)

    mrsim.simulate_data()

    input_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'
    fpath_output = input_fpath_prefix + 'mr_static_simulation.h5'

    output_file = Path(fpath_output)
    if not output_file.is_file():
        mrsim.write_simulation_results(fpath_output)

    return 1


def test_motion_mr_simulation(rec=False, verb=False, throw=True):

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/cDynamicSimulation/'

    fpath_xml = input_fpath_prefix + 'Segmentations/XCAT_TissueParameters_XML.xml'
    fpath_template_rawdata = input_fpath_prefix + 'TemplateData/MR/CV_nav_cart_64Cube_1Echo.h5'

    #
    mrsim, rawdata, labels = prepare_test_simulation(fpath_template_rawdata, fpath_xml)
    
    mrsim.set_acquisition_template_data(rawdata)

    csm = pMR.CoilSensitivityData()
    csm.calculate(rawdata)
    mrsim.set_csm(csm)

    SNR = 5
    SNR_label = 13

    mrsim.set_snr(SNR)
    mrsim.set_snr_label(SNR_label)

    # 
    num_resp_states = 4
    resp_motion = pDS.MRMotionDynamic(num_resp_states)
    resp_motion.set_cyclicality(False)

    # generate artificial motion signal
    Nt = 100
    t0_s = 0
    tmax_s = 1200
    time_points = np.linspace(t0_s, tmax_s, Nt)

    resp_frequency_Hz = 0.2
    resp_curve = 0.5 * ( 1 + np.sin( 2*np.pi*resp_frequency_Hz*time_points))

    resp_motion.set_dynamic_signal(time_points, resp_curve)

    #
    inhale_dvf = prep_displacement_field(labels)
    identity_trafo = prep_displacement_field(labels)
    identity_trafo.fill(0)
    
    resp_motion.add_displacement_field(identity_trafo)
    resp_motion.add_displacement_field(inhale_dvf)
    

    #
    resp_motion.set_mr_acquisitions(rawdata)

    # 
    mrsim.add_motion_dynamic(resp_motion)
    mrsim.simulate_data()

    #   
    input_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'
    fpath_output = input_fpath_prefix + 'mr_motion_simulation.h5'

    output_file = Path(fpath_output)
    if not output_file.is_file():
        mrsim.write_simulation_results(fpath_output)

    return 1

def test_main(rec=False, verb=False, throw=True):
    
    num_tests = 0
    # num_tests += test_static_mr_simulation(rec, verb, throw)
    num_tests += test_motion_mr_simulation(rec, verb, throw)

    return False, num_tests

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
