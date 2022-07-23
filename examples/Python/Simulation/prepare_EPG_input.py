import sys 
sys.path.append("/home/sirfuser/devel/buildVM/sources/SIRF/src/xDynamicSimulation/pDynamicSimulation/")

import numpy as np

import TissueParameterList as TPL

root_path = '/media/sf_CCPPETMR/TestData/Input/xDynamicSimulation/pDynamicSimulation/' 
xml_path = root_path + 'Cube128/XCAT_TissueParameters_XML.xml'
fpath_output = root_path + 'Fingerprints/XCAT_tissue_parameter_list.npz'


# parse file
tpl = TPL.TissueParameterList()
tpl.parse_xml(xml_path)
tpl.print_contents()

# prepare EPG input array
all_mr_params = tpl.mr_as_array()
mr_params = all_mr_params[:,1:]
mr_params[:,0] /= 100
mr_params[:,-1] = 0

mr_params_unique,idx_inverse = np.unique(mr_params, axis=0, return_inverse=True)

np.savez(fpath_output, mr_params_full=all_mr_params, mr_parameters=mr_params_unique, unique_idx_inverse=idx_inverse)
