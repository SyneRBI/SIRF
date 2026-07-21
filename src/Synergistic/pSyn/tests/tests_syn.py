import pytest
from sirf.Reg import NiftiImageData3D
from sirf.SIRF import ImageData
from sirf.Utilities import examples_data_path, existing_filepath


@pytest.mark.parametrize("engine", ["Gadgetron", "STIR"])
def test_nifti_image_data_3d(engine):
    filename = {"Gadgetron": "SIRF_recon.h5", "STIR": "dicom_as_nifti.nii"}[engine]
    data_path = examples_data_path('MR')
    raw_data_file = existing_filepath(f"{data_path}/zenodo", filename)
    image = ImageData()
    image.read(raw_data_file, engine, 1)
    nifti_image = NiftiImageData3D(image)
    assert image == nifti_image
