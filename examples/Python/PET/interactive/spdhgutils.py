#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
###
# Utility functions for
# Demonstration of advanced PET reconstruction with CCP PET-MR Software
# and randomized algorithms
#
# This demonstration shows how to use the stochastic primal-dual hybrid
# gradient algorithm (SPDHG) for regularized PET reconstruction.
#
#    [CERS2018] A. Chambolle, M. J. Ehrhardt, P. Richtarik and C.-B. Schoenlieb,
#    *Stochastic Primal-Dual Hybrid Gradient Algorithm with Arbitrary Sampling
#    and Imaging Applications*. SIAM Journal on Optimization, 28(4), 2783–2808
#    (2018) http://doi.org/10.1007/s10851-010-0251-1
#
#    [E+2017] M. J. Ehrhardt, P. J. Markiewicz, P. Richtarik, J. Schott,
#    A. Chambolle and C.-B. Schoenlieb, *Faster PET reconstruction with a
#    stochastic primal-dual hybrid gradient method*. Wavelets and Sparsity XVII,
#    58 (2017) http://doi.org/10.1117/12.2272946.
#
#    [EMS2018] M. J. Ehrhardt, P. J. Markiewicz and C.-B. Schoenlieb, *Faster
#    PET Reconstruction with Non-Smooth Priors by Randomization and
#    Preconditioning*. (2018) ArXiv: http://arxiv.org/abs/1808.07150
#
# This demo is a 'script', i.e. intended to be run step by step in a
# Python IDE such as spyder. It is organised in 'cells'. spyder displays these
# cells nicely and allows you to run each cell on its own.
#
# First version: 8th of September 2018
# Author: Matthias J Ehrhardt, Edoardo Pasca
#
## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
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
import numpy
import pSTIR as pet
from ccpi.plugins.regularisers import FGP_TV, TGV, LLT_ROF, Diff4th

class OperatorInd(object):

    def __init__(self, op, subset_num, num_subsets):
        self.__op__ = op
        self.__subset_num__ = subset_num
        self.__num_subsets__ = num_subsets

        # FIXME: this may be known to STIR
        x = op.img_templ.copy()
        x.fill(1)
        y = self.forward_sirf(x).as_array()
        self.ind = numpy.nonzero(y.flatten())[0]

    def forward_sirf(self, x):
        return self.__op__.forward(x, subset_num=self.__subset_num__,
                                   num_subsets=self.__num_subsets__)

    def __call__(self, x):
        y = self.__op__.forward(x, subset_num=self.__subset_num__,
                                num_subsets=self.__num_subsets__)
        return self.sirf2sub(y)

    def direct(self, x):
        y = self.__op__.direct(x, subset_num=self.__subset_num__,
                                num_subsets=self.__num_subsets__)
        return self.sirf2sub(y)

    def forward(self, x):
        return self(x)

    def adjoint(self, x):
        x = self.sub2sirf(x)
        return self.__op__.adjoint(x, subset_num=self.__subset_num__,
                                    num_subsets=self.__num_subsets__)

    def allocate_direct(self, x=None):
        y = self.__op__.img_templ.copy()
        if x is not None:
            y.fill(x)
        return y

    def allocate_adjoint(self, x=None):
        y = pet.AcquisitionData(self.__op__.acq_templ)
        if x is not None:
            y.fill(x)
        return self.sirf2sub(y)

    def sub2sirf(self, x):
        y = pet.AcquisitionData(self.__op__.acq_templ)
        y.fill(0)
        y_array = y.as_array().flatten()
        y_array[self.ind] = x
        y.fill(y_array)
        return y

    def sirf2sub(self, x):
        #y = numpy.zeros(len(self.ind))
        return x.as_array().flatten()[self.ind]

#class SubsetOperator():
#
#    def __init__(self, op, nsubsets):
#        self.__op__ = op
#        self.__nsubsets__ = nsubsets
#
#    def __len__(self):
#        return self.__nsubsets__
#
#    def __call__(self, x):
#        return self.__op__.forward(x)
#
#    def adjoint(self, x):
#        return self.__op__.backward(x)
#
#    def __ind__(self, i):
#        return OperatorInd(self.__op__, i, len(self))


def SubsetOperator(op, nsubsets):
    return [OperatorInd(op, ind, nsubsets) for ind in range(nsubsets)]

# the FGP_TV will output a CCPi DataContainer not a SIRF one, so
# we will need to wrap it in something compatible
class FGP_TV_SIRF(FGP_TV):
    '''SIRF wrapper for CCPi CIL plugin FGP_TV'''
    def prox(self, x, sigma):
        print("calling FGP")
        out = super(FGP_TV_SIRF, self).prox(x, sigma)
        print("done")
        y = x.copy()
        y.fill(out.as_array())
        return y

class TGV_TV_SIRF(TGV):
    '''SIRF wrapper for CCPi CIL plugin TGV'''
    def prox(self, x, sigma):
        print("calling FGP")
        out = super(TGV_TV_SIRF, self).prox(x, sigma)
        print("done")
        y = x.copy()
        y.fill(out.as_array())
        return y

class ROF_LLT_SIRF(LLT_ROF):
    '''SIRF wrapper for CCPi CIL plugin LLT_ROF'''
    def prox(self, x, sigma):
        print("calling ROF_LLT")
        out = super(ROF_LLT_SIRF, self).prox(x, sigma)
        print("done")
        y = x.copy()
        y.fill(out.as_array())
        return y

class Diff4th_SIRF(Diff4th):
    '''SIRF wrapper for CCPi CIL plugin Diff4th'''
    def prox(self, x, sigma):
        print("calling Diff4th")
        out = super(Diff4th_SIRF, self).prox(x, sigma)
        print("done")
        y = x.copy()
        y.fill(out.as_array())
        return y
def PowerMethodNonsquare(op, numiters, x0=None):
    # Initialise random
    # Jakob's
    #inputsize = op.size()[1]
    #x0 = ImageContainer(numpy.random.randn(*inputsize)
    # Edo's
    #vg = ImageGeometry(voxel_num_x=inputsize[0],
    #                   voxel_num_y=inputsize[1],
    #                   voxel_num_z=inputsize[2])
    #
    #x0 = ImageData(geometry = vg, dimension_labels=['vertical','horizontal_y','horizontal_x'])
    #print (x0)
    #x0.fill(numpy.random.randn(*x0.shape))
    print ("call powermethod")

    if x0 is None:
        x0 = op.create_image_data()

    #s = numpy.zeros(numiters)
    # Loop
    for it in numpy.arange(numiters):
        x1 = op.adjoint(op.direct(x0))
        x1norm = numpy.sqrt(x1.dot(x1))
        #print ("x0 **********" ,x0)
        #print ("x1 **********" ,x1)
        #s[it] = x1.dot(x0) / x0.dot(x0)
        x0 = (1.0/x1norm)*x1
    return numpy.sqrt(x1norm)
