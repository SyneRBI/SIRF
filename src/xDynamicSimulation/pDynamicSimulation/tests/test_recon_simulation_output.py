# -*- coding: utf-8 -*-
"""sirf.Gadgetron test.
v{version}

Constructor of MR image data from MR acquisition data

Usage:
  test_imagedata_constructor [--help | options]

Options:
  -v, --verbose  report each test status
{author}

{licence}
"""

import numpy as np 

from sirf.Gadgetron import *
from sirf.Utilities import runner, __license__
__version__ = "0.2.3"
__author__ = "Johannes Mayer"


def test_recon_output_simulate_dynamics(record=False, verb=False, throw=True):

    print("Running a reconstruction of simulated MR data")

    prefix_data_path = "/media/sf_CCPPETMR/TestData/Output/xDynamicSimulation/"
    input_data_path = prefix_data_path + "cDynamicSimulation/"
    
    rawdata = AcquisitionData(input_data_path + '/output_test_test_simulate_dynamics.h5')
    print("We have {} acquisitions.".format(rawdata.number()))
    rawdata = preprocess_acquisition_data(rawdata)
    rawdata.sort()
    
    imgdata = ImageData()
    imgdata.from_acquisition_data(rawdata)
        
    acq_model = AcquisitionModel()
    acq_model.set_up(rawdata, imgdata)

    csms = CoilSensitivityData()
    csms.calculate(rawdata)
    acq_model.set_coil_sensitivity_maps(csms)

    recon = acq_model.inverse(rawdata)

    dicom_value_range = (2**16-1)
    img_content = recon.as_array()
    img_content =  dicom_value_range * (img_content - np.amin(img_content[:]) ) / ( np.amax(img_content[:]) - np.amin(img_content[:]) )

    recon.fill(img_content)

    output_data_path = prefix_data_path + "pDynamicSimulation/"

    recon.write(output_data_path + "output_recon_simulate_dynamics.dcm")

    test_failed = False
    return test_failed, 1


def test_main(record=False, verb=False, throw=True):
    test_recon_output_simulate_dynamics(record, verb, throw)


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
