import argparse
import numpy
import os
import pylab
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
AcquisitionData demo
''')
args = parser.parse_args()

def main():

    # direct all information printing to a file
    info_printer = printerTo('test_ad_info.txt', INFO_CHANNEL)
    # write all warnings in a file
    warn_printer = printerTo('test_ad_warn.txt', WARNING_CHANNEL)

    # define acquisition data
    ad = AcquisitionData('my_raw_data.hs')
    ada = ad.as_array()
    print(ada.shape)
    ac = ad.clone()
    aca = ac.as_array()
    print(aca.shape)
    aca -= ada
    print(numpy.amax(aca))
    ac = ad.get_empty_copy(1.0)
    aca = ac.as_array()
    print(numpy.amax(aca))

try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
