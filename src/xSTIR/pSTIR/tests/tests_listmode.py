# -*- coding: utf-8 -*-
"""sirf.STIR tests
v{version}

Usage:
  tests_listmode [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import sirf.STIR as pet
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Richard Brown"


def test_main(rec=False, verb=False, throw=True):
    msg_red = pet.MessageRedirector()

    data_path = pet.examples_data_path('PET')
    raw_data_file = pet.existing_filepath(data_path, 'mMR/list.l.hdr')

    lm2sino = pet.ListmodeToSinograms()
    lm2sino.set_input(raw_data_file)
    
    prompt_rate_threshold = 73036.
    known_time = 22.
    time_at_which_prompt_rate_exceeds_threshold = \
        lm2sino.get_time_at_which_prompt_rate_exceeds_threshold(prompt_rate_threshold)

    if abs(time_at_which_prompt_rate_exceeds_threshold-known_time) > 1.e-4:
        raise AssertionError("ListmodeToSinograms::get_time_at_which_prompt_rate_exceeds_threshold failed")


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
