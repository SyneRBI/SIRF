'''MR test sets.

Usage:
  test_all [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

record = args['--record']
verbose = args['--verbose']

import sys

failed = 0
ntests = 0

#print('\n\n--- fully sampled reconstruction test set:')
#import fully_sampled
#failed += fully_sampled.main()

#print('\n\n--- undersampled reconstruction test set:')
#import undersampled
#failed += undersampled.main()

print('\n\n--- test set 1:')
import test1
failed1, ntests1 = test1.main()
failed += failed1
ntests += ntests1

print('\n\n--- test set 2:')
import test2
failed2, ntests2 = test2.main()
failed += failed2
ntests += ntests2

print('\n\n--- test set 3:')
import test3
failed3, ntests3 = test3.main()
failed += failed3
ntests += ntests3

if record:
    print('%d measurements recorded' % ntests)
    sys.exit(0)

if failed == 0:
    print('all tests passed')
    sys.exit(0)

print('%d of the tests failed' % failed)
sys.exit(failed)
