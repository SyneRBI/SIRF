#!/usr/bin/env python

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
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

import glob
import os
import sys

for i in glob.glob('*.py'):
    narg = len(sys.argv)
    if narg > 1 and i.find('listmode') >= 0:
        continue
    if os.path.abspath(__file__) == os.path.abspath(i):
        continue
    print('\n=== %s\n' % i)
    args = ''
    for a in range(1, narg):
        args += ' ' + sys.argv[a]
    exe = '"' + sys.executable + '" ' + i + args
    err = os.system(exe)
    if err:
        raise RuntimeError(i + ' failed')
