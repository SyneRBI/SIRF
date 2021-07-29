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

import sirf.DynamicSimulation as pDS
import sirf.Gadgetron as pMR
import sirf.Reg as pReg

from sirf.Utilities import __license__, runner

__version__ = "3.1.0"
__author__ = "Johannes Mayer"


def test_main(rec=False, verb=False, throw=True):

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/cDynamicSimulation/'

    fpath_xml = input_fpath_prefix + 'Segmentations/XCAT_TissueParameters_XML.xml'
    fpath_template_rawdata = input_fpath_prefix + 'TemplateData/MR/CV_nav_cart_64Cube_1Echo.h5'

    rawdata = pMR.AcquisitionData(fpath_template_rawdata)
    rawdata = pMR.preprocess_acquisition_data(rawdata)

    empty_img = pMR.ImageData()
    empty_img.from_acquisition_data (rawdata)

    labels = pReg.NiftiImageData3D(empty_img)
    label_array = labels.as_array()
    max_label = 50
    label_array = np.random.randint(1, max_label,size = label_array.shape)
    labels.fill(label_array)

    mrsim = pDS.DynamicSimulation(labels, fpath_xml)
    mrsim.set_acquisition_template_data(rawdata)

    csm = pMR.CoilSensitivityData()
    csm.calculate(rawdata)

    mrsim.set_csm(csm)

    SNR = 15
    SNR_label = 1

    mrsim.set_snr(SNR)
    mrsim.set_snr_label(SNR_label)

    mrsim.simulate_data()

    input_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'
    fpath_output = input_fpath_prefix + 'mr_static_simulation.h5'
    mrsim.write_simulation_results(fpath_output)

    return False, 1

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
