'''
Demo for CSMs computation from coil images by xGadgetron and ismrmrd-python-tools
'''

import argparse
import math
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'
DATA_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/examples/'
IPT_PATH = os.environ.get('SRC_PATH') + '/ismrmrd-python-tools/ismrmrdtools'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)
sys.path.append(IPT_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Demo for CSMs computation from coil images by xGadgetron and ismrmrd-python-tools
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()

try:
 
    from ismrmrdtools import coils, show

    input_data = MR_Acquisitions(DATA_PATH + args.filename)
    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    print('sorting acquisitions...')
    processed_data.sort()

    CIs = MR_CoilImages()
    print('computing coil images...')
    CIs.calculate(processed_data)

    CSMs = MR_CoilSensitivityMaps()
    print('computing sensitivity maps...')
    CSMs.calculate(CIs, method = 'SRSS(niter = 10)')

    nz = CSMs.number()
    print('%d slices' % nz)

    while HAVE_PYLAB:
        print('---\n Enter the slice number to view it.')
        print(' A value outside the range [1 : %d] will stop this loop.'% nz)
        s = str(input('slice: '))
        if len(s) < 1:
            break
        z = int(s) - 1
        if z < 0 or z >= nz:
            break

        coil_images = numpy.squeeze(CIs.as_ndarray(z))

        print('computing sensitivity maps (Inati)...')
        (csm_inati, rho) = coils.calculate_csm_inati_iter\
                           (coil_images, niter = 10)

        csm_srss = numpy.squeeze(CSMs.as_ndarray(z))
        nc, ny, nx = csm_srss.shape
        images = numpy.ndarray((2, ny, nx), dtype = numpy.float64)
        images[0, :, :] = abs(numpy.sum(csm_srss * coil_images, axis = 0))
        images[1, :, :] = abs(numpy.sum(csm_inati * coil_images, axis = 0))
        maxv = numpy.amax(images)
        show.imshow(images, tile_shape = (1,2), scale = (0, maxv),\
            titles = ['Square Root of the Sum of Squares', 'Inati'])

        csms = numpy.ndarray((2*nc, ny, nx), dtype = numpy.float64)
        csms[ 0 :   nc, :, :] = abs(csm_srss)
        csms[nc : 2*nc, :, :] = abs(csm_inati)
        t = ['SRSS']
        t.extend(['' for x in range(nc - 1)])
        t.extend(['Inati'])
        t.extend(['' for x in range(nc - 1)])
        show.imshow(csms, tile_shape = (4,4), scale = (0,1), titles = t)

        print('computing sensitivity maps (Walsh)...')
        (csm_walsh, rho) = coils.calculate_csm_walsh(coil_images)
        csms[nc : 2*nc, :, :] = abs(csm_walsh)
        t[8] = 'Walsh'
        show.imshow(csms, tile_shape = (4,4), scale = (0,1), titles = t)

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
