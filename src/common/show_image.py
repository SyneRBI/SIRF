'''
Usage:
  show_image.py <image_file>
'''

import sys

narg = len(sys.argv)
if narg < 2:
    print('\nUsage:\n')
    print('python show_image.py <image_file>')
image_file = sys.argv[1]

i = image_file.find('.h5')
if i > 0:
    from sirf.Gadgetron import ImageData
    image_data = ImageData(image_file)
    image_data.show(title='Image data', cmap=None)
else:
    from sirf.STIR import ImageData
    image_data = ImageData(image_file)
    image_data.show(title='Image data')
