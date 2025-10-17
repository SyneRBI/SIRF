# -*- coding: utf-8 -*-
"""sirf.Gadgetron Views Test.
v{version}

DataContainer algebra tests

Usage:
  data_views.py

{author}

{licence}
"""
from sirf.SIRF import norm, dot, copyto
from sirf.Gadgetron import AcquisitionData, \
    AcquisitionDataView, ImageDataView, FullySampledReconstructor
from sirf.Utilities import examples_data_path, error
__version__ = "3.1.0"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def main():
    data_path = examples_data_path('MR')
    AcquisitionData.set_storage_scheme('memory')
    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    x = processed_data
    y = x + 0
    z = x + 0

    x_view = AcquisitionDataView(x)
    y_view = AcquisitionDataView(y)
    z_view = AcquisitionDataView(z)

    y = 2*x
    norm_y = y.norm()
    copyto(y_view, x_view)
    y_view *= 2
    norm_y_view = norm(y_view)
    print(norm_y, norm_y_view)
    z = x + y
    norm_z = z.norm()
    copyto(z_view, x_view)
    z_view += y_view
    norm_z_view = norm(z_view)
    print(norm_z, norm_z_view)
    s = norm(x)
    t = norm(x_view)
    print(s, t)
    s = x.sum()
    t = x_view.sum()
    print(s, t)
    s = x.dot(y)
    t = dot(x_view, y_view)
    print(s, t)

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    x = complex_images
    y = x + 0
    z = x + 0

    x_view = ImageDataView(x)
    y_view = ImageDataView(y)
    z_view = ImageDataView(z)

    y = 2*x
    norm_y = y.norm()
    copyto(y_view, x_view)
    y_view *= 2
    norm_y_view = norm(y_view)
    print(norm_y, norm_y_view)
    z = x + y
    norm_z = z.norm()
    copyto(z_view, x_view)
    z_view += y_view
    norm_z_view = norm(z_view)
    print(norm_z, norm_z_view)
    s = norm(x)
    t = norm(x_view)
    print(s, t)
    s = x.sum()
    t = x_view.sum()
    print(s, t)
    s = x.dot(y)
    t = dot(x_view, y_view)
    print(s, t)

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)

