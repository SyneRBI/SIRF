import sys

from sirf.SIRF import ImageData
from sirf.Reg import NiftiImageData3D

na = len(sys.argv)
if na < 3:
    sys.exit(1)

image = ImageData()
image.read(sys.argv[1], sys.argv[2], 1)
nifti_image = NiftiImageData3D(image)
sys.exit(image != nifti_image)