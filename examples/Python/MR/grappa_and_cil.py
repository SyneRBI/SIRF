'''
GRAPPA reconstruction with an iterative algorithm from CIL: illustrates
the use of AcquisitionModel in CIL optimisation 

Usage:
  grappa_and_cil.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2019 University College London.
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

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import sirf
from sirf.Utilities import existing_filepath
from sirf.Utilities import error
from sirf.Utilities import show_3D_array
from sirf.Gadgetron import petmr_data_path
from sirf.Gadgetron import AcquisitionData, ImageData
from sirf.Gadgetron import AcquisitionModel
from sirf.Gadgetron import AcquisitionDataProcessor
from sirf.Gadgetron import CartesianGRAPPAReconstructor
from sirf.Gadgetron import CoilSensitivityData


from ccpi.optimisation.funcs import Norm2sq
from ccpi.optimisation.funcs import ZeroFun
from ccpi.optimisation.algs import FISTA, FBPD, CGLS
#from ccpi.optimisation.ops import PowerMethodNonsquare
from ccpi.plugins.regularisers import FGP_TV#, TGV, LLT_ROF, Diff4th
from ccpi.framework import DataContainer as cilDataContainer

import numpy
import time

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('mr')
class Algorithm(object):
    
    def __init__(self):
        self.iteration = 0
        self.stop_cryterion = 'max_iter'
        self.__max_iteration = 0
        self.__loss = []
        self.memopt = False
        self.timing = []
    def set_up(self, *args, **kwargs):
        raise NotImplementedError()
    def update(self):
        raise NotImplementedError()
    
    def should_stop(self):
        '''stopping cryterion'''
        raise NotImplementedError()
    
    def __iter__(self):
        return self
    def next(self):
        '''python2 backwards compatibility'''
        return self.__next__()
    def __next__(self):
        if self.should_stop():
            raise StopIteration()
        else:
            time0 = time.time()
            self.update()
            self.timing.append( time.time() - time0 )
            self.iteration += 1
    def get_output(self):
        '''Returns the solution found'''
        return self.x
    def get_current_loss(self):
        '''Returns the current value of the loss function'''
        return self.__loss[-1]
    @property
    def loss(self):
        return self.__loss
    @property
    def max_iteration(self):
        return self.__max_iteration
    @max_iteration.setter
    def max_iteration(self, value):
        assert isinstance(value, int)
        self.__max_iteration = value
    
class GradientDescent(Algorithm):
    '''Implementation of a simple Gradient Descent algorithm
    '''
    
    def __init__(self, **kwargs):
        '''initialisation can be done at creation time if all 
        proper variables are passed or later with set_up'''
        super(GradientDescent, self).__init__()
        self.x = None
        self.rate = 0
        self.objective_function = None
        self.regulariser = None
        args = ['x_init', 'objective_function', 'rate']
        for k,v in kwargs.items():
            if k in args:
                args.pop(args.index(k))
        if len(args) == 0:
            return self.set_up(x_init=kwargs['x_init'],
                               objective_function=kwargs['objective_function'],
                               rate=kwargs['rate'])
    
    def should_stop(self):
        '''stopping cryterion, currently only based on number of iterations'''
        return self.iteration >= self.max_iteration
    
    def set_up(self, x_init, objective_function, rate):
        '''initialisation of the algorithm'''
        self.x = x_init.copy()
        if self.memopt:
            self.x_update = x_init.copy()
        self.objective_function = objective_function
        self.rate = rate
        self.loss.append(objective_function(x_init))
        
    def update(self):
        '''Single iteration'''
        if self.memopt:
            self.objective_function.gradient(self.x, out=self.x_update)
            self.x_update *= -self.rate
            self.x += self.x_update
        else:
            self.x += -self.rate * self.objective_function.grad(self.x)
            
        self.loss.append(self.objective_function(self.x))
        


class FISTAAlg(Algorithm):
    '''Fast Iterative Shrinkage-Thresholding Algorithm
    
    Beck, A. and Teboulle, M., 2009. A fast iterative shrinkage-thresholding 
    algorithm for linear inverse problems. 
    SIAM journal on imaging sciences,2(1), pp.183-202.
    
    Parameters:
      x_init: initial guess
      f: data fidelity
      g: regularizer
      h:
      opt: additional algorithm 
    '''
    
    def __init__(self, **kwargs):
        '''initialisation can be done at creation time if all 
        proper variables are passed or later with set_up'''
        super(FISTAAlg, self).__init__()
        self.f = None
        self.g = None
        self.invL = None
        self.t_old = 1
        args = ['x_init', 'f', 'g', 'opt']
        for k,v in kwargs.items():
            if k in args:
                args.pop(args.index(k))
        if len(args) == 0:
            return self.set_up(x_init=kwargs['x_init'],
                               f=kwargs['f'],
                               g=kwargs['g'],
                               opt=kwargs['opt'])
    
    def set_up(self, x_init, f=None, g=None, opt=None):
        
        # default inputs
        if f   is None: 
            self.f = ZeroFun()
        else:
            self.f = f
        if g   is None:
            g = ZeroFun()
        else:
            self.g = g
        
        # algorithmic parameters
        if opt is None: 
            opt = {'tol': 1e-4, 'iter': 1000, 'memopt':False}
        
        self.max_iteration = opt['iter'] if 'iter' in opt.keys() else 1000
        self.tol = opt['tol'] if 'tol' in opt.keys() else 1e-4
        memopt = opt['memopt'] if 'memopt' in opt.keys() else False
        self.memopt = memopt
            
        # initialization
        if memopt:
            self.y = x_init.clone()
            self.x_old = x_init.clone()
            self.x = x_init.clone()
            self.u = x_init.clone()
        else:
            self.x_old = x_init.copy()
            self.y = x_init.copy()
        
        #timing = numpy.zeros(max_iter)
        #criter = numpy.zeros(max_iter)
        
    
        self.invL = 1/f.L
        
        self.t_old = 1
        
    def should_stop(self):
        '''stopping cryterion, currently only based on number of iterations'''
        return self.iteration >= self.max_iteration
    
    def update(self):
    # algorithm loop
    #for it in range(0, max_iter):
    
        if self.memopt:
            # u = y - invL*f.grad(y)
            # store the result in x_old
            self.f.gradient(self.y, out=self.u)
            self.u.__imul__( -self.invL )
            self.u.__iadd__( self.y )
            # x = g.prox(u,invL)
            self.g.proximal(self.u, self.invL, out=x)
            
            self.t = 0.5*(1 + numpy.sqrt(1 + 4*(self.t_old**2)))
            
            # y = x + (t_old-1)/t*(x-x_old)
            self.x.subtract(self.x_old, out=self.y)
            self.y.__imul__ ((self.t_old-1)/self.t)
            self.y.__iadd__( self.x )
            
            self.x_old.fill(self.x)
            self.t_old = self.t
            
            
        else:
            u = self.y - self.invL*self.f.grad(self.y)
            
            self.x = self.g.prox(u,self.invL)
            
            self.t = 0.5*(1 + numpy.sqrt(1 + 4*(self.t_old**2)))
            
            self.y = self.x + (self.t_old-1)/self.t*(self.x-self.x_old)
            
            self.x_old = self.x.copy()
            self.t_old = self.t
        
        
        self.loss.append( self.f(self.x) + self.g(self.x) )
        
class cilPluginToSIRFFactory(object):
    '''Factory to create SIRF wrappers for CCPi CIL plugins'''
    @staticmethod
    def getInstance(thetype, **kwargs):
        '''Returns an instance of a CCPi CIL plugin wrapped to work on SIRF DataContainers'''
        obj = thetype(**kwargs)
        orig_prox = obj.prox
        obj.prox = cilPluginToSIRFFactory.prox(orig_prox, 
                                               obj.__class__.__name__)
        return obj
    @staticmethod
    def prox(method, classname):
        def wrapped(x, sigma):
            '''Wrapped method'''
            print("calling ", classname)
            if isinstance(x, sirf.Gadgetron.ImageData):
                # if the data is MR => complex we operate the regulariser
                # only on the real part
                X = x.as_array()
                out = method(cilDataContainer(X.real), sigma)
                X.real = out.as_array()[:]
                y = x.copy()
                y.fill(X)
            else:
                out = method(x, sigma)
                y = x.copy()
                y.fill(out)
            print("done")
                
            return y
        return wrapped

    
from ccpi.optimisation.funcs import Function
class SumFunction(Function):
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2
    def __call__(self, x):
        return self.f1(x) + self.f2(x)
    def gradient(self, x , out = None):
        return self.f1.gradient(x) + self.f2.gradient(x)
    def grad(self, x):
        return self.gradient(x, None)
    
def PowerMethodNonsquare(op,numiters , x0=None):
    # Initialise random
    
    if x0 is None:
        #x0 = op.create_image_data()
        x0 = op.allocate_direct()
        x0.fill(numpy.random.randn(*x0.shape))
    
    s = numpy.zeros(numiters)
    # Loop
    for it in numpy.arange(numiters):
        x1 = op.adjoint(op.direct(x0))
        #x1norm = numpy.sqrt((x1*x1).sum())
        x1norm = x1.norm()
        #print ("x0 **********" ,x0)
        #print ("x1 **********" ,x1)
        s[it] = (x1*x0.conjugate()).sum() / (x0*x0.conjugate()).sum()
        x0 = (1.0/x1norm)*x1
    return numpy.sqrt(s[-1]), numpy.sqrt(s), x0


if True:

    # locate the input data file
    input_file = existing_filepath(data_path, data_file)

    # acquisition data will be read from an HDF file input_data
    acq_data = AcquisitionData(input_file)

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    # the following function runs a GadgetChain data processor
    # preprocessed_data = preprocess_acquisition_data(acq_data)
    gadget_chain = ['NoiseAdjustGadget', \
                    'AsymmetricEchoAdjustROGadget', \
                    'RemoveROOversamplingGadget']
    ap = AcquisitionDataProcessor( gadget_chain )
    ap.set_input( acq_data )
    ap.process()
    preprocessed_data = ap.get_output()


    # perform reconstruction
    recon = CartesianGRAPPAReconstructor()
    recon.set_input( preprocessed_data )
    recon.compute_gfactors( False )
    print('---\n reconstructing...')
    recon.process()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    image_data = recon.get_output()
    
    image_array = image_data.as_array()
    title = 'Reconstructed image data (magnitude)'
    show_3D_array(abs(image_array), suptitle = title, label = 'slice', \
                  xlabel = 'samples', ylabel = 'readouts', show = False)
    
    # compute coil sensitivity maps
    csms = CoilSensitivityData()
    print('---\n sorting acquisition data...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate( preprocessed_data )

    # create acquisition model based on the acquisition parameters
    # stored in preprocessed_data and image parameters stored in
    # image_data
    acq_model = AcquisitionModel( preprocessed_data, image_data )
    acq_model.set_coil_sensitivity_maps( csms )

    # use the acquisition model (forward projection) to simulate acquisition data
    simulated_data = acq_model.forward( image_data )

    # create a Norm2square objective function
    # c || Ax - y ||^2
    norm2sq = Norm2sq( A = acq_model , b = simulated_data , c = 1)
    # create a random initialisation image by shuffling the real
    # image data. 
    x_init = image_data.copy()
    x = x_init.as_array().flatten()
    numpy.random.shuffle(x)
    x = numpy.reshape(x, x_init.as_array().shape)
    x_init.fill(x)
    del x
    show_3D_array(x_init.as_array().real)
    
    # test if <Ax0,y0> = <y0, A^Tx0>
    y0 = simulated_data.copy()
    x = y0.as_array().flatten()
    numpy.random.shuffle(x)
    x = numpy.reshape(x, y0.as_array().shape)
    y0.fill(x)
    del x     
    x0 = x_init
    fx0 = acq_model.direct(x0)
    by0 = acq_model.adjoint(y0)
    a = fx0.dot(y0)
    b = by0.dot(x0)
    numpy.testing.assert_almost_equal(abs((a-b)/a), 0, decimal=4)
    
    
    # calculate Lipschitz constant
    lipschitz = PowerMethodNonsquare( acq_model , numiters = 10 , x0 = x_init) [0] 
    norm2sq.L = lipschitz  
    print ("Lipschitz " , norm2sq.L)
    
#%%
    # create a Gradient Descent algorithm which minimises norm2sq
    gd = GradientDescent(x_init=x_init, 
               objective_function=norm2sq, rate=lipschitz/3.)

    gd.max_iteration = 20
    pixval = []
    gadgval = image_data.as_array()[0][46][160]
    import matplotlib.pyplot as plt
    for i,el in enumerate(gd):
    #for i in range(10):
    #    gd.update()
        pixval.append( gd.get_output().as_array()[0][46][160])
        if i%1 == 0:
            print ("\rIteration {} Loss: {} pix {}".format(gd.iteration, 
               gd.get_current_loss(), pixval[-1]/gadgval))
    
    fig = plt.figure()
    plt.plot([i/gadgval for i in gd.loss])
    plt.show()
    
    #show_3D_array(gd.get_output().as_array().real, suptitle='Gradient Descent (magnitude)')

#%%
    # USE FISTA with Regularisation   
    no_regulariser = ZeroFun()
    # create a regulariser with the Factory
    regulariser = cilPluginToSIRFFactory.getInstance(FGP_TV, 
                                           lambdaReg=1.,
                                           iterationsTV=300,
                                           tolerance=1e-5,
                                           methodTV=0,
                                           nonnegativity=0,
                                           printing=0,
                                           device='cpu')
    options = {'tol': 1e-4, 'iter': 10, 'memopt':False}

    norm2sq.L = lipschitz*3.
    # create a FISTA algorithm instance
    fista = FISTAAlg(x_init=x_init, f=norm2sq, g=regulariser, opt=options)
    fpixval = []
    #%%
    # run FISTA
    for i,el in enumerate(fista):
        fpixval.append( fista.get_output().as_array()[0][46][160])
        if i%1 == 0:
            print ("\rFISTA Iteration {} Loss: {} pix {}".format(fista.iteration, 
               fista.get_current_loss(), fpixval[-1]/gadgval))
#%%
    fista_noreg = FISTAAlg(x_init=x_init, f=norm2sq, g=no_regulariser, opt=options)
    fpixval = []
    #%%
    # run FISTA
    for i,el in enumerate(fista_noreg):
        fpixval.append( fista_noreg.get_output().as_array()[0][46][160])
        if i%1 == 0:
            print ("\rFISTA Iteration {} Loss: {} pix {}".format(fista_noreg.iteration, 
               fista_noreg.get_current_loss(), fpixval[-1]/gadgval))

    #refined_image_array = fista.get_output().as_array()

    # show reconstructed and refined images
    
    title = 'FISTA image data (magnitude)'
    #show_3D_array(abs(refined_image_array), suptitle = title, label = 'slice', \
    #              xlabel = 'samples', ylabel = 'readouts')
    fig = plt.figure()
    #plt.imshow(gd.get_output().as_array()[0].real)
    #plt.show()
    plt.plot([i/gadgval for i in fista.loss])
    plt.show()
    
    #%%
    # plot the results
    fig = plt.figure()
    ax1 = plt.subplot(1,4,1)
    plt.imshow(abs(image_data.as_array()[0]))
    ax1.set_title('Initial Data')
    ax2 = plt.subplot(1,4,2)
    plt.imshow(abs(gd.get_output().as_array()[0]))
    ax2.set_title('Gradient Descent')
    ax2 = plt.subplot(1,4,3)
    plt.imshow(abs(fista_noreg.get_output().as_array()[0]))
    ax2.set_title('FISTA no reg')
    ax2 = plt.subplot(1,4,4)
    plt.imshow(abs(fista.get_output().as_array()[0]))
    ax2.set_title('FISTA + FGP_TV')
    plt.show()
#%%
#    options['iter'] = 20
#    out = FBPD(x_init, operator=acq_model, constraint=None, data_fidelity=norm2sq,\
#         regulariser=None, opt=options)
#    # x, it, timing, criter
#    title = 'FBPD'
#    show_3D_array(abs(out[0].as_array()), suptitle = title, label = 'slice', \
#                  xlabel = 'samples', ylabel = 'readouts')
    #%%
#    options['iter'] = 50
#    # x_init, operator , data , opt=None):
#    outcgls = CGLS(x_init, operator=acq_model, data=simulated_data, opt=options)
#    # x, it, timing, criter
#    title = 'CGLS'
#    show_3D_array(abs(outcgls[0].as_array()), suptitle = title, label = 'slice', \
#                  xlabel = 'samples', ylabel = 'readouts')
#try:
#    main()
#    print('done')
#
#except error as err:
#    # display error information
#    print('??? %s' % err.value)
