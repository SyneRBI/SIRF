# -*- coding: utf-8 -*-
"""sirf.STIR.AcquisitionModel tests
v{version}

Usage:
  test_six_adjoint.py [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import sirf.STIR as pet
from sirf.Utilities import is_operator_adjoint, runner, __license__
__version__ = "0.2.3"
__author__ = "Ander Biguri"

def test_main(rec=False, verb=False, throw=True):

    pet.MessageRedirector()

    for scheme in ("file", "memory"):
        pet.AcquisitionData.set_storage_scheme(scheme)

        original_verb = pet.get_verbosity()
        pet.set_verbosity(False)
        # create an acq_model that is explicitly a RayTracingMatrix
        am = pet.AcquisitionModelUsingRayTracingMatrix()
        # load sample data
        data_path = pet.examples_data_path('PET')
        raw_data_file = pet.existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
        ad = pet.AcquisitionData(raw_data_file)
        # create sample image
        image = pet.ImageData()
        image.initialise(dim=(31, 111, 111), vsize=(2.25, 2.25, 2.25))
        # set up Acquisition Model
        am.set_up(ad,image)
        # test for adjointnesss
        if not is_operator_adjoint(am, verbose = verb):
          raise AssertionError('AcquisitionModelUsingRayTracingMatrix is not adjoint')

    # Reset original verbose-ness
    pet.set_verbosity(original_verb)
    return 0, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
