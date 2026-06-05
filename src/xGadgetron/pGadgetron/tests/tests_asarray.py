"""sirf.Gadgetron.ImageData.asarray tests.

v{version}

Usage:
  tests_asarray.py [--help | options]

Options:
  -r, --record   not used
  -v, --verbose  report each test status

{author}

{licence}
"""
import numpy

from sirf.SIRF import ContiguousError
from sirf.Utilities import runner, RE_PYEXT, examples_data_path, pTest

__version__ = '0.1.0'
__author__ = "Evgueni Ovtchinnikov"


def test_main(rec=False, verb=False, throw=True, no_ret_val=True, \
        data_path=examples_data_path('MR'), engine='Gadgetron'):
    # import engine module
    import importlib

    mr = importlib.import_module('sirf.' + engine)
    mr.AcquisitionData.set_storage_scheme('memory')

    # process command-line options
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    print('\n-- testing discontiguous image:')
    import os.path
    acq_data = mr.AcquisitionData(os.path.join(data_path, '..', 'MR', 'simulated_MR_2D_cartesian.h5'))
    preprocessed_data = mr.preprocess_acquisition_data(acq_data)
    recon = mr.FullySampledReconstructor()
    recon.set_input(preprocessed_data)
    recon.process()
    img_data = recon.get_output()
    try:
        test.ntest += 1
        img_data.asarray(copy=False)
    except ContiguousError:
        pass
    else:
        test.failed = True
        print('expected ContiguousError not raised')

    try:
        test.ntest += 1
        diff = img_data.asarray() - img_data.as_array()
        test.failed = numpy.linalg.norm(diff) > 0
    except Exception as e:
        test.failed = True
        print(e)

    numpy.testing.assert_equal(test.failed, 0)
    if no_ret_val:
        return
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
