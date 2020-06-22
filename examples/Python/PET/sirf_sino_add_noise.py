# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2020 University College London.
#
# This is software developed for the Collaborative Computational
# Project in Positron Emission Tomography and Magnetic Resonance imaging
# (http://www.ccppetmr.ac.uk/).
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from sirf.Utilities import error
import sirf.STIR as pet
import numpy as np
import argparse

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __contains__(self, item):
        return self.__eq__(item)

    def __iter__(self):
        yield self

    def __str__(self):
        return '[{0},{1}]'.format(self.start, self.end)

parser = argparse.ArgumentParser(description='Add noise to a sinogram.')
parser.add_argument('percentage', type=float, help='Percentage of counts', metavar='percentage', choices=[Range(0.0, 100.0)])
parser.add_argument('sino_in', type=str, help='Input sinogram')
parser.add_argument('sino_out', type=str, nargs='?', help='Output sinogram prefix')
args = parser.parse_args()

def main():
    sino = pet.AcquisitionData(args.sino_in)
    sino = add_poisson_noise(args.percentage/100., sino)

    if args.sino_out:
        f_out = args.sino_out
    else:
        f_out = args.sino_in + "-" + args.percentage.replace('.','_') + '-percent'
    sino.write(f_out)


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
