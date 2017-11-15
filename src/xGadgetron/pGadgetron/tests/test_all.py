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

import sys

failed = 0

print('\n\n--- fully sampled reconstruction test set:')
import fully_sampled
failed += fully_sampled.main()

print('\n\n--- undersampled reconstruction test set:')
import undersampled
failed += fully_sampled.main()

print('\n\n--- test set 3:')
import test3
failed3, ntests3 = test3.main()
failed += failed3

if failed == 0:
    print('all tests passed')
    sys.exit(0)
else:
    print('%d of the tests failed' % failed)
    sys.exit(failed)
