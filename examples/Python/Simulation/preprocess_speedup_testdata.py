
import numpy as np 
from pathlib import Path
import nibabel as nib



num_slices = 10
Nz = 128
slab_start = Nz//2 -num_slices//2
slab_end = Nz//2 + num_slices//2


def preprocess_mvfs(rootpath, folder_pattern, prefix_output):

    fpath_input = Path(rootpath + folder_pattern)
    print("looking in {} ".format(fpath_input))
    list_files = sorted(fpath_input.glob('mvf_*'))

    for f in list_files:
        print("loading {}".format(f))
        
        mvf = nib.load(str(f))
        print("The input motoinfield has size {}".format(mvf.shape))

        mvf = mvf.slicer[:,:,slab_start:slab_end,:]
        print("The output motoinfield has size {}".format(mvf.shape))
        
        fname_output = str(prefix_output + folder_pattern + f.name)
        print("Storing to {}".format(fname_output))
        nib.save(mvf, fname_output)

root_path = '/media/sf_CCPPETMR/TestData/Input/xDynamicSimulation/pDynamicSimulation/'

fpath_input = root_path + 'Cube128/'
fpath_output = root_path + 'Slab128/'

foldername_resp = 'mvf_resp/'
foldername_card = 'mvf_card/'

preprocess_mvfs(fpath_input, foldername_resp, fpath_output)
preprocess_mvfs(fpath_input, foldername_card, fpath_output)


fname_segmentation = fpath_input + 'label_volume.nii'
seg = nib.load(fname_segmentation)

print("The input segmentation has size {}".format(seg.shape))

seg = seg.slicer[:,:,slab_start:slab_end]
nib.save(seg, fpath_output + 'label_volume.nii')