'''
Usage:
  show_image.py <image_file>
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2020 Rutherford Appleton Laboratory STFC
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


import sys


def show(image_file):
    i = image_file.find('.h5')
    if i > 0:
        from sirf.Gadgetron import ImageData
        image_data = ImageData(image_file)
        image_data.show(title='Image data', cmap=None)
    else:
        from sirf.STIR import ImageData
        image_data = ImageData(image_file)
        image_data.show(title='Image data')


if __name__ == "__main__":
    narg = len(sys.argv)
    if narg < 2:
        print('\nUsage:\n')
        print('python show_image.py <image_file>')
    else:
        show(sys.argv[1])