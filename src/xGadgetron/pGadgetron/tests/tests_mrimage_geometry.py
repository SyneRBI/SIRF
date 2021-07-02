
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

import sirf.Reg as sreg
from sirf.Gadgetron import *
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Johannes Mayer"


def test_2d_slice_stack(rec=False, verb=False, throw=True):

    print("Running the img geometry test")
    data_path = examples_data_path('MR')
    rawdata = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    rawdata = preprocess_acquisition_data(rawdata)
    rawdata.sort()
    
    imgdata = ImageData()
    imgdata.from_acquisition_data(rawdata)
    
    print("The image positions are {:.2f}".format(imgdata.get_ISMRMRD_info('position')))

    print("The rawdata directions are: {:.2f}".format(rawdata.get_ISMRMRD_info('read_dir')[0]))
    print("The rawdata directions are: {:.2f}".format(rawdata.get_ISMRMRD_info('phase_dir')[0]))
    print("The rawdata directions are: {:.2f}".format(rawdata.get_ISMRMRD_info('slice_dir')[0]))

    print("The image data directions are: {:.2f}".format(imgdata.get_ISMRMRD_info('read_dir')[0]))
    print("The image data directions are: {:.2f}".format(imgdata.get_ISMRMRD_info('phase_dir')[0]))
    print("The image data directions are: {:.2f}".format(imgdata.get_ISMRMRD_info('slice_dir')[0]))

    recon = FullySampledReconstructor()
    recon.set_input(rawdata)
    recon.process()
    img_data = recon.get_output()
    
    img_data = img_data.abs()
    img_data.write( data_path+ "tmp_mrgeometry.dcm")
    
    nii_img = sreg.NiftiImageData(img_data)
    img_data.write( data_path+ "tmp_mrgeometry.nii")
    
    
    test_failed = False
    return test_failed, 1


if __name__ == "__main__":
    runner(test_2d_slice_stack, __doc__, __version__, __author__)
