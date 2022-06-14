''' 
Object-Oriented wrap for the cGadgetron-to-Python interface pygadgetron.py
'''

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
# Copyright 2018 - 2020 University College London
# Copyright 2018 - 2022 Physikalisch-Technische Bundesanstalt (PTB)
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
    '''
    Class to perform dynamic MRI simulations
    '''
    def __init__(self, tissue_labels, fname_xml):

        self.handle = None
        assert_validity(tissue_labels, pReg.NiftiImageData3D)

        self.handle = pysim.cDS_MRDynamicSimulation(tissue_labels.handle, fname_xml)
        check_status(self.handle)

    def set_template_data(self, ad):
        '''
        Method to set both contrast and acquisition template data, cf. methods set_acquisition_template_data and set_contrast_template_data. 
        '''
        self.set_acquisition_template_data(ad)
        self.set_contrast_template_data(ad)

    def set_acquisition_template_data(self, ad):
        '''
        Setting MR rawdata to define
        - a destination geometry into which the signal-filled segmentation is resampled.
        - defining the encoding that is performed, i.e. trajectory and # of readout and phase-encoding points etc.
        '''
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setAcquisitionTemplateData(self.handle, ad.handle)


    def set_contrast_template_data(self, ad):
        '''
        Setter to define MR rawdata to fix sequence parameters used in signal computation for steady state sequences.
        '''
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setContrastTemplateData(self.handle, ad.handle)

    
    def simulate_data(self):
        '''
        Method to start simulation after it is set up.
        '''
        try_calling( pysim.cDS_simulateData(self.handle) )


    def set_csm(self, csm):
        '''
        Setter for the employed coil sensitiviy profiles. They have to match the geometry of the acquisition template data.
        '''
        assert_validity(csm, pMR.CoilSensitivityData)
        try_calling( pysim.cDS_setCoilmaps(self.handle, csm.handle))
    
    def set_snr(self, SNR):
        '''
        Setting the signal to noise ratio which is added to the MR rawdata as complex Gaussian noise
        '''
        SNR = np.array(SNR).astype(np.float32)
        pysim.cDS_setSNR(self.handle, SNR.ctypes.data)


    def set_snr_label(self, SNR_label):
        '''
        Method to define a label in the segmentation corresponding to a tissue type used to compute the signal from on which the SNR is based.
        '''
        pysim.cDS_setNoiseLabel(self.handle, SNR_label)


    def set_offset_trafo(self, trafo):
        '''
        Method to add an affine transform which is applied to the 3D segmentation prior to resampling into template geometry
        '''
        assert_validity(trafo, pReg.AffineTransformation)
        try_calling(pysim.cDS_setOffsetTransformation(self.handle, trafo.handle))

    def write_simulation_results(self, fpath_output):
        pysim.cDS_writeSimulationResults(self.handle, fpath_output) 

    def save_motion_ground_truth(self):
        pysim.cDS_saveMotionGroundTruth(self.handle)


    def save_parametermap_ground_truth(self, filename_prefix):
        '''
        Method storing the parameters T1, T2, spin density, chemical shift and the segmentation resampled into the destination geometry as nifti images.
        '''
        pysim.cDS_saveParameterMapsGroundTruth(self.handle, filename_prefix)
        return (filename_prefix + "_T1_ms.nii", filename_prefix + "_T2_ms.nii", filename_prefix + "_spindensity.nii", filename_prefix + "_labels.nii")
    
    '''
    Method to include an independent motion mode in the simulation.
    '''
    def add_motion_dynamic(self, motiondyn):
        try_calling(pysim.cDS_addMRMotionDynamic(self.handle, motiondyn.handle)) 

    '''
    Method to include contrast changes during the data acquisition.
    '''
    def add_contrast_dynamic(self, contdyn):
        try_calling(pysim.cDS_addMRContrastDynamic(self.handle, contdyn.handle))

    '''
    Method to activate the simulation of a pre-computed magnetisation. This will overwrite other contrast dynamics.
    '''
    def add_external_contrast_dynamic(self, contrastdyn):
        try_calling(pysim.cDS_addExternalContrastDynamic(self.handle, contrastdyn.handle))
        

class Dynamic(object):
    '''
    Class to define properties shared by dynamic processes during data acquisition simulation that can be added to the simulation class.
    '''
    def __init__(self, num_states):
        self.num_states = num_states

    
    def set_dynamic_signal(self, time_points_seconds, signal_points):
        '''
        Method to set the temporal evolution of the simulation 
        '''
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

    def get_idx_corr_sizes(self, ad):

        assert_validity(ad, pMR.AcquisitionData)
        sizes = np.zeros(shape=(self.num_states,), dtype=np.int32)
        sizes = np.ascontiguousarray(sizes)
        pysim.cDS_getIdxCorrSizes(self.handle, ad.handle, sizes.ctypes.data)
        
        return sizes

    def get_idx_corr(self, ad):
        '''
        Method to get the index of the phase encoding lines that were simulated in the same dynamic state.
        '''
        idx_sizes = self.get_idx_corr_sizes(ad)
        idx_corr = []
        for ibin in range(len(idx_sizes)):
            num_entries=idx_sizes[ibin]
            idx = np.zeros(shape=(num_entries,), dtype=np.int32)
            idx = np.ascontiguousarray(idx)
            pysim.cDS_getIdxCorr(self.handle, ibin, idx.ctypes.data)
            idx_corr.append(idx)
        
        return idx_corr


class MRMotionDynamic(Dynamic):
    '''
    Class to define motion that is executed during the simulted data acquisition.
    '''
    def __init__(self, num_states):
        '''
        num_states: # of motion states that are simulated. This is NOT constrained to be equal to the # of motion vector fields 
                    that are supplied to the motion dynamic.
        '''
        self.num_states = num_states
        self.handle = None
        self.handle = pysim.cDS_MRMotionDynamic(num_states)
        check_status(self.handle)
   

    def add_displacement_field(self, dvf):
        '''
        Method to add a displacement field, each added displacement vector field  defines a motion state. The order in which the 
        motion fields are supplied is important where the first corresponds to surrogate signal = 0 and the last to surrogate signal = 1.
        '''
        assert_validity(dvf, pReg.NiftiImageData3DDisplacement)
        pysim.cDS_addMRDisplacementField(self.handle, dvf.handle)


    def set_cyclicality(self, is_cyclic):
        '''
        Method defining whether surrogate signal = 0 corresponds to the same motion state as surrogate signal = 1.
        '''
        pysim.cDS_setCyclicality(self.handle, is_cyclic)

    

    def set_groundtruth_folder_prefix(self, prefix_existing_path):
        '''
        Method defining the path to which the motion dynamic stores its ground truth motion fields.
        '''
        pysim.cDS_setMRGroundTruthFolderName(self.handle, prefix_existing_path)

    

class TissueParameter():
    '''
    Class carrying a tissue parameter associated with a simulation.
    The information is filled based on the XML file detailing the .
        
    This is primarily used to define tissue parameters corresponding to surrogate signal = 0 and surrogate signal = 1
    in an MRContrastDynamic.
    '''
    def __init__(self, sim, label):
        '''
        sim: MRDynamicSimulation object
        label: label in segmentation corresponding to the tissue.
        '''
        self.handle = None

        assert_validity(sim, MRDynamicSimulation)
        
        self.handle = pysim.cDS_getTissueParameter(sim.handle, label)
        check_status(self.handle)

    def set_T1_value(self, T1_ms):
        '''
        Method to modify the T1 value of the tissue parameters.
        '''
        pysim.cDS_setT1Value(self.handle,T1_ms)

    def set_spin_density(self, spin_density):
        '''
        Method to modify the spin density of the tissue parameters.
        '''
        pysim.cDS_setSpinDensity(self.handle, spin_density)



class MRContrastDynamic(Dynamic):
    '''
    Class to modify the tissue parameters over time during the simulation of MR data acquisition.
    This is primarily used to modify T1 contrast agent inflow.
    '''
    def __init__(self, num_states):
        '''
        num_states: # of simualted contrast states across the time of data acquisition.
        '''
        self.num_states = num_states
        self.handle = None
        self.handle = pysim.cDS_MRContrastDynamic(num_states)
        check_status(self.handle)

    def add_dynamic_label(self, label):
        '''
        Method to add segmentation labels that will behave identially during the simulation of contrast variation.
        '''
        pysim.cDS_addDynamicLabel(self.handle, label)

    def set_parameter_extremes(self, tissue_0, tissue_1):
        '''
        tissue_0: tissue parameter corresponding to surrogate signal = 0.
        tissue_1: tissue parameter corresponding to surrogate signal = 1.
        '''
        assert_validity(tissue_0, TissueParameter)
        assert_validity(tissue_1, TissueParameter)

        pysim.cDS_setMRParameterExtremes(self.handle, tissue_0.handle, tissue_1.handle)

class ExternalMRContrastDynamic(Dynamic):
    '''
    Class to carry information on pre-computed magnetisation that is simulated.
    '''
    def __init__(self):
        self.handle = None
        self.handle = pysim.cDS_ExternalMRContrastDynamic()
        check_status(self.handle)
    
    def add_external_signal(self, external_signal):
        
        num_sig_pts = external_signal.get_num_signal_points()
        num_labels = external_signal.get_num_labels()
        print("Adding external MR signal for {} labels at {} time points".format(num_labels, num_sig_pts))

        for i in range(num_sig_pts):
            ptr_labels, ptr_sig = external_signal.get_signal_pointers(i)
            pysim.cDS_appendExternalTissueSignal(self.handle, num_labels, ptr_labels, ptr_sig)

# helper class to set up external MR signal
class ExternalMRSignal():
    '''
    Auxiliary class for storing information on pre-computed magnetisation. Each label in a segmentation must have an 
    associated magnetisation signal that is available for each simulated time point.
    '''
    def __init__(self, labels, signals):
        '''
        labels: comprehensive set of labels to which the magnetisation in signals
        signals: shape (labels.size, #time points) where # time points = acquisition_template.number()
        '''
        assert type(signals) is np.ndarray, "Please pass a numpy ndarray of signals. You provided {}".format(type(signals))
        assert np.iscomplexobj(signals), "Please pass a 64-bit complex numpy array"

        (num_tissue_signals, num_time_pts) = signals.shape
        self.number_of_labels = labels.size

        assert self.number_of_labels == num_tissue_signals, "Please pass a signal array of shape (#labels, #time points)"

        self.tissue_signals = np.array([self.TissueSignal() for _ in range(num_time_pts)])
        self.set_up_signal(labels, signals)

    def get_num_signal_points(self):
        return self.tissue_signals.size

    def get_num_labels(self):
        return self.number_of_labels

    def set_up_signal(self, labels, signals):
        
        if signals.dtype is not np.complex64:
            the_signals = signals.astype(np.complex64)
        else:
            the_signals = signals

        if labels.dtype is not np.int32:
            the_labels = labels.astype(np.int32)
        else:
            the_labels = labels

        num_time_pts = the_signals.shape[1]

        for t in range(num_time_pts):
            signal_at_time = the_signals[:,t]
            self.tissue_signals[t] = self.TissueSignal(the_labels, signal_at_time)
    '''
    Auxiliary method to access arrays in C-level.
    '''
    def get_signal_pointers(self, time_index):

        return self.tissue_signals[time_index].get_signal_pointers()
    '''
    Auxiliary class to access arrays in C-level.
    '''
    class TissueSignal():

        def __init__(self, labels=None, signals=None):
            
            self.time_ms = 0.0
            self.labels = np.ascontiguousarray(labels)
            self.signals = np.ascontiguousarray(signals)
        
        def get_signal_pointers(self):
            return( self.labels.ctypes.data, self.signals.ctypes.data)