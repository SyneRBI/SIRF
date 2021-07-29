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
# Created on Tue Nov 21 10:17:28 2017

import sirf.DynamicSimulation as pDS
import sirf.Gadgetron as pMR
import sirf.Reg as pReg

from sirf.Utilities import __license__, runner

__version__ = "3.1.0"
__author__ = "Johannes Mayer"


def test_main(rec=False, verb=False, throw=True):

    fpath_prefix = '/media/sf_CCPPETMR/TestData/Input/xDynamicSimulation/cDynamicSimulation/'
    
    fpath_xml = fpath_prefix + 'Segmentations/XCAT_TissueParameters_XML.xml'

    fpath_template_rawdata = fpath_prefix + 'TemplateData/MR/CV_nav_cart_64Cube_1Echo.h5'

    rawdata = pMR.AcquisitionData(fpath_template_rawdata)
    rawdata = pMR.preprocess_acquisition_data(rawdata)

    empty_img = pMR.ImageData()
    empty_img.from_acquisition_data (rawdata)

    labels = pReg.NiftiImageData3D(empty_img)
    label_array = labels.as_array()
    label_array[:] = 1
    labels.fill(label_array)

    sim = pDS.DynamicSimulation(labels, fpath_xml)
    sim.set_acquisition_template_data(rawdata)

    
    return False, 1
    

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
