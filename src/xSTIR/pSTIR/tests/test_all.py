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

print('\n\n--- test set 1:')
import test1
failed1, ntest1 = test1.main()

print('\n\n--- test set 2:')
import test2
failed2, ntest2 = test2.main(verb = True)

print('\n\n--- test set 3:')
import test3
failed3, ntest3 = test3.main(verb = True)

print('\n\n--- test set 4:')
import test4
failed4, ntest4 = test4.main(verb = True)

failed = failed1 + failed2 + failed3 + failed4
ntest = ntest1 + ntest2 + ntest3 + ntest4

if failed == 0:
    print('all %d tests passed' % ntest)
    sys.exit(0)
else:
    print('%d of %d tests failed' % (failed, ntest))
    sys.exit(failed)
