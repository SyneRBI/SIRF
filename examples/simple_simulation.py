'''
Upper-level demo that illustrates the computation of coil sensitivity maps,
applying projection from the image space into acquisition space and back
defined by the aquisition model, and images and acquisitions algebra.
'''

import argparse
import os
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level demo that illustrates the computation of coil sensitivity maps,
applying projection from the image space into acquisition space and back
defined by the aquisition model, and images and acquisitions algebra.
''')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    print('---\n acquisition data norm: %e' % input_data.norm())

    print('---\n processing acquisitions...')
    processed_data = input_data.process(['RemoveROOversamplingGadget'])

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction
    recon = SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % complex_images.norm())

    csms = CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(processed_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in complex_images
    am = AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce 'acquisitions'
    acqs = am.forward(complex_images)

    print('---\n reconstructed images forward projection norm %e' % acqs.norm())

    # get data as a Python ndarray
    acqs_data = acqs.as_array();

    # TODO display a slice etc

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
