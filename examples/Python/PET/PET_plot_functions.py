'''
Some functions for plotting used in the examples
'''

## CCP SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2020-2021 University College London
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

import sirf.STIR as PET
import numpy as np
import matplotlib.pyplot as plt

def plot_sinogram_profile(prompts, randoms=None, scatter=None, sumaxis=(0,1), select=0):
    """
    Plot a profile through sirf.STIR.AcquisitionData

    sumaxis: axes to sum over (passed to numpy.sum(..., axis))
    select: element to select after summing
    """
    if isinstance(prompts, str):
        prompts = PET.AcquisitionData(prompts)
    if isinstance(randoms, str):
        randoms = PET.AcquisitionData(randoms)
    if isinstance(scatter, str):
        scatter = PET.AcquisitionData(scatter)
    # we will average over all sinograms to reduce noise
    plt.figure()
    ax = plt.subplot(111)
    plt.plot(np.sum(prompts.as_array(), axis=sumaxis)[select,:], label='prompts')
    if randoms is None:
        if scatter is not None:
            plt.plot(np.sum(scatter.as_array(), axis=sumaxis)[select,:], label='scatter')
    else:
        randoms_as_array = randoms.as_array()
        plt.plot(np.sum(randoms_as_array, axis=sumaxis)[select,:], label='randoms')
        if scatter is not None:
            plt.plot(np.sum(scatter.as_array() + randoms_as_array, axis=sumaxis)[select,:], label='randoms+scatter')
    ax.legend()
    plt.show()

