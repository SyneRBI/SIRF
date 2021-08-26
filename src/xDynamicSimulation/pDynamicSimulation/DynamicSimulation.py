''' 
Object-Oriented wrap for the cGadgetron-to-Python interface pygadgetron.py
'''

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
# Copyright 2018 - 2020 University College London
# Copyright 2018 - 2021 Physikalisch-Technische Bundesanstalt (PTB)
#
# This is software developed for the Collaborative Computational
# Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
# (http://www.ccpsynerbi.ac.uk/).
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import abc
import numpy as np
from numbers import Number, Complex, Integral


try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False

import sys

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

from deprecation import deprecated

from sirf.Utilities import check_status, assert_validity, try_calling, error

import sirf
from sirf import SIRF
from sirf.SIRF import DataContainer

import sirf.pysimulation as pysim
import sirf.pyiutilities as pyiutil
# import sirf.pygadgetron as pygadgetron
# import sirf.pysirf as pysirf

import sirf.Reg as pReg
import sirf.Gadgetron as pMR



def format_arrays_for_setters(data):

    if not isinstance(data, np.ndarray):
        raise error('Wrong input format.' + \
            ' Should be numpy.ndarray. Got {}'.format(type(data)))

    if data.dtype != np.float32:
            the_data = data.astype(np.float32)

    convert = not the_data.flags['C_CONTIGUOUS']
    if convert:
        the_data = np.ascontiguousarray(the_data)

    return the_data


class MRDynamicSimulation(object):

    def __init__(self, tissue_labels, fname_xml):

        self.handle = None
        assert_validity(tissue_labels, pReg.NiftiImageData3D)

        self.handle = pysim.cDS_MRDynamicSimulation(tissue_labels.handle, fname_xml)
        check_status(self.handle)

    def set_acquisition_template_data(self, ad):
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setAcquisitionTemplateData(self.handle, ad.handle)

    def set_contrast_template_data(self, ad):
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setContrastTemplateData(self.handle, ad.handle)

    def simulate_data(self):
        pysim.cDS_simulateData(self.handle)

    def set_csm(self, csm):
        assert_validity(csm, pMR.CoilSensitivityData)
        pysim.cDS_setCoilmaps(self.handle, csm.handle)

    def set_snr(self, SNR):
        pysim.cDS_setSNR(self.handle, SNR)

    def set_snr_label(self, SNR_label):
        pysim.cDS_setNoiseLabel(self.handle, SNR_label)

    def write_simulation_results(self, fpath_output):
        pysim.cDS_writeSimulationResults(self.handle, fpath_output); 

    def save_motion_ground_truth(self):
        pysim.cDS_saveMotionGroundTruth(self.handle)

    def add_motion_dynamic(self, motiondyn):
        pysim.cDS_addMRMotionDynamic(self.handle, motiondyn.handle); 

class Dynamic(object):

    def set_dynamic_signal(self, time_points_seconds, signal_points):
      
        num_signal_points = time_points_seconds.size
        if num_signal_points != signal_points.size:
            AssertionError("The signal and time point arrays do not have identical size.")

        the_time = format_arrays_for_setters(time_points_seconds)
        the_signal = format_arrays_for_setters(signal_points)
        
        pysim.cDS_setDynamicSignal(self.handle, 
                                    the_time.ctypes.data,\
                                    the_signal.ctypes.data,\
                                    num_signal_points)

    def add_displacement_field(self, dvf):
        assert_validity(dvf, pReg.NiftiImageData3DDisplacement)
        pysim.cDS_addMRDisplacementField(self.handle, dvf.handle)

    def set_cyclicality(self, is_cyclic):
        pysim.cDS_setCyclicality(self.handle, is_cyclic)

 
class MRMotionDynamic(Dynamic):

    def __init__(self, num_states):
        self.handle = None
        self.handle = pysim.cDS_MRMotionDynamic(num_states)
        check_status(self.handle)
   
    def add_displacement_field(self, dvf):
        assert_validity(dvf, pReg.NiftiImageData3DDisplacement)
        pysim.cDS_addMRDisplacementField(self.handle, dvf.handle)

    def set_cyclicality(self, is_cyclic):
        pysim.cDS_setCyclicality(self.handle, is_cyclic)

    def set_groundtruth_folder_prefix(self, prefix_existing_path):
        pysim.cDS_setMRGroundTruthFolderName(self.handle, prefix_existing_path)