
# -*- coding: utf-8 -*-
"""sirf.Gadgetron test.
v{version}

Constructor of MR image data from MR acquisition data

Usage:
  test_imagedata_constructor [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status
{author}

{licence}
"""

import numpy as np 

from sirf.Gadgetron import *
from sirf.Utilities import runner, __license__
__version__ = "0.2.3"
__author__ = "Johannes Mayer"

def test_MR_DICOM_writer(rec=False, verb=False, throw=True):

    print("Running the img geometry test")
    data_path = examples_data_path('MR')
    rawdata = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    rawdata = preprocess_acquisition_data(rawdata)
    rawdata.sort()
   
    recon = FullySampledReconstructor()
    recon.set_input(rawdata)
    recon.process()
    img_data = recon.get_output()
    
    print("The image positions are {}".format(img_data.get_ISMRMRD_info('position')))

    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('read_dir')[0]))
    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('phase_dir')[0]))
    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('slice_dir')[0]))

    print("The image data directions are: {}".format(img_data.get_ISMRMRD_info('read_dir')[0]))
    print("The image data directions are: {}".format(img_data.get_ISMRMRD_info('phase_dir')[0]))
    print("The image data directions are: {}".format(img_data.get_ISMRMRD_info('slice_dir')[0]))

    img_data = img_data.abs()

    img_content = img_data.as_array()
    dicom_value_range = (2**16-1)
    img_content =  dicom_value_range * (img_content - np.amin(img_content[:]) ) / ( np.amax(img_content[:]) - np.amin(img_content[:]) )
    img_data.fill(img_content)

    img_data.write( data_path+ "/output_test_MR_DICOM_writer.dcm")
    
    test_failed = False
    return test_failed, 1

if __name__ == "__main__":
    runner(test_MR_DICOM_writer, __doc__, __version__, __author__)
