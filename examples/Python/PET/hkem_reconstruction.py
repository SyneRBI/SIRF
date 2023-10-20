'''Hybrid Kernelized Expectation Maximization reconstruction demo.
Implements a Kernelized Ordered Subsets (OS) version of the One Step Late
algorithm (OSL) from Green et al for Maximum a Posteriori (MAP) maximisation.

Usage:
  hkem_reconstruction [--help | options]

Options:
  -f <file>, --file=<file>     raw data file [default: my_forward_projection.hs]
  -a <file>, --anim=<file>     anatomical image file [default: test_image_PM_QP_6.hv]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -s <subs>, --subs=<subs>     number of subsets [default: 12]
  -i <iter>, --subiter=<iter>  number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  --non-interactive            do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2023 National Physical Laboratory
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London.
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

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from cgitb import small
from matplotlib.pyplot import title
import numpy as np

from sirf.Utilities import error, examples_data_path, existing_filepath

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)


# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
if args['--anim'] is not None:
    ai_file = existing_filepath(data_path, args['--anim'])
else:
    ai_file = None
show_plot = not args['--non-interactive']

def divide(numerator, denominator, small_num):
    """ division like in STIR 
    """
    small_value = np.max(numerator)*small_num
    if small_value<=0:
        small_value=0
    X,Y,Z= numerator.shape
    for i in range(Z):
        for j in range(Y):
            for k in range(X):
                if numerator[k][j][i]<=small_value and denominator[k][j][i]<=small_value:
                   numerator[k][j][i]=0
                else:
                    numerator[k][j][i]/=denominator[k][j][i]
    return numerator

def divide_sino(numerator, denominator, small_num):
    """ division like in STIR 
    """
    small_value = np.max(numerator)*small_num
    if small_value<=0:
        small_value=0
    X,Y,Z= numerator.shape
    for i in range(Z):
        for j in range(Y):
            for k in range(X):
                if numerator[k][j][i]<=small_value and denominator[k][j][i]<=small_value:
                   numerator[k][j][i]=0
                else:
                    numerator[k][j][i]/=denominator[k][j][i]
    return numerator


# Define a function that does something with an image. This function
# provides a simplistic example of user's involvement in the reconstruction
def image_data_processor(image_array, im_num):
    """ Process/display an image

        image is not modified in this simplistic example - but might have been
    """
    if not show_plot:
        return image_array
    # display the current estimate of the image at z = 20
    pylab.figure(im_num)
    pylab.title('image estimate %d' % im_num)
    pylab.imshow(image_array[20,:,:])
    print('You may need to close Figure %d window to continue' % im_num)
    return image_array


def main():

    # direct all engine's information and warnings printing to files
    msg_red = pet.MessageRedirector('info.txt', 'warn.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from this file
    # (TODO: a link to raw data formats document to be given here)
    print('raw data: %s' % raw_data_file)
    acq_data = pet.AcquisitionData(raw_data_file)

    # read anatomical image
    anatomical_image = pet.ImageData(ai_file)
    if show_plot:
        anatomical_image.show(title='Image used as anatomical prior')
    image_array = anatomical_image.as_array()
    image_array[image_array < 0] = 0
    anatomical_image.fill(image_array)

    # create initial image estimate
    image = anatomical_image.get_uniform_copy()
    current_alpha1 = anatomical_image.get_uniform_copy(1)
    current_alpha2 = anatomical_image.get_uniform_copy(1)
    iterative_kernel_info1 = current_alpha1.get_uniform_copy(1)
    iterative_kernel_info2 = current_alpha2.get_uniform_copy(1)
    image_update1 = current_alpha1.get_uniform_copy(1)
    image_update2 = current_alpha1.get_uniform_copy(1)
    if show_plot:
        image.show(title='Image used as initial guess')

    # set up acquisition model
    acq_model.set_up(acq_data, image)
    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)

    # select Kernelized Ordered Subsets Maximum A-Posteriori One Step Late
    # as the reconstruction algorithm
    recon = pet.KOSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(1)
    recon.set_num_subiterations(num_subiterations)
    recon.set_input(acq_data)
    recon.set_anatomical_prior(anatomical_image)
    recon.set_num_neighbours(5)
    recon.set_num_non_zero_features(1)
    recon.set_sigma_m(2.0)
    recon.set_sigma_p(3.0)
    recon.set_sigma_dm(5.0)
    recon.set_sigma_dp(5.0)
    recon.set_only_2D(True)
    recon.set_hybrid(False)

    # set up the reconstructor based on a sample image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    sensistivity = acq_model.backward(acq_data.get_uniform_copy(1))
    mult_update = image.get_uniform_copy(1)
    diff_im=image.get_uniform_copy()
    # now reconstruct an image using the KOSMAPOSL update function
    recon.set_up(current_alpha1)
    # set the initial image estimate
    sens1=obj_fun.get_subset_sensitivity(0)
    # alternatively:
    #sens1=recon.get_subset_sensitivity()

    Ksensitivity1 = recon.compute_kernelised_image(sens1, iterative_kernel_info1)
    recon.set_current_estimate(current_alpha1)
    for subiteration in range(1):
        print('\n------------- sub-iteration %d' % subiteration)
        # perform one KOSMAPOSL sub-iteration
        recon.update_current_estimate()
        current_alpha1 = recon.get_current_estimate()
        image_update1 = recon.compute_kernelised_image(current_alpha1, iterative_kernel_info1)   
        iterative_kernel_info1 = current_alpha1  
        recon.set_current_estimate(current_alpha1)  

# Now let's create a KOSMAPOSL algorithm in python
    recon2 = pet.KOSMAPOSLReconstructor()
    recon2.set_objective_function(obj_fun)
    recon2.set_num_subsets(1)
    recon2.set_num_subiterations(num_subiterations)
    recon2.set_input(acq_data)
    recon2.set_anatomical_prior(anatomical_image)
    recon2.set_num_neighbours(5)
    recon2.set_num_non_zero_features(1)
    recon2.set_sigma_m(2.0)
    recon2.set_sigma_p(3.0)
    recon2.set_sigma_dm(5.0)
    recon2.set_sigma_dp(5.0)
    recon2.set_only_2D(True)
    recon2.set_hybrid(False)
    recon2.set_up(current_alpha2)
    recon2.set_current_estimate(current_alpha2)
# in order to see the reconstructed image evolution
    # open up the user's access to the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for subiteration in range(1):#num_subiterations
        print('\n------------- sub-iteration %d' % subiteration)
        # perform one KOSMAPOSL sub-iteration
        Kalpha = recon2.compute_kernelised_image(current_alpha2, iterative_kernel_info2)
        Ksensitivity = recon2.compute_kernelised_image(sensistivity, iterative_kernel_info2)
        gradient_plus_sensitivity = acq_model.backward(acq_data/acq_model.forward(Kalpha))# G = B (Y/F* K*a)
        Kgradient_plus_sensitivity = recon2.compute_kernelised_image(gradient_plus_sensitivity, iterative_kernel_info2) #K*G
        mult_update_array = divide(Kgradient_plus_sensitivity.as_array(),Ksensitivity.as_array(),0.000001)#K*G/S
        mult_update.fill(mult_update_array)
        current_alpha2 *=mult_update #a(n+1)=a(n) * K*G/S
        image_update2=recon.compute_kernelised_image(current_alpha2, iterative_kernel_info2) #Lambda(n+1)=K(n)a(n+1)
        iterative_kernel_info2=current_alpha2 #K(n+1)<=K(n)
    
    mean2=np.mean(image_update2.as_array())
    mean1=np.mean(image_update1.as_array())
    diff_im.fill(abs(image_update2.as_array()/mean2-image_update1.as_array()/mean1))
    max_diff=np.max(diff_im.as_array())
    mean_diff=np.mean(diff_im.as_array())
    max_im=np.max(image_update1.as_array())
    rel_diff=(max_diff)/(max_im)*100
    print('Max perc difference : %e' % (rel_diff))
    if rel_diff < 1:
        print('Max difference is less then 1% so is probably OK')
    else:
        print('Max difference is higher then 1% something could be wrong')

    # forward projection of the reconstructed image simulates the
    # acquisition of data by the scanner
    print('projecting...')
    simulated_data = acq_model.forward(image)
    # compute the reconstruction residual
    diff = simulated_data * (acq_data.norm()/simulated_data.norm()) - acq_data
    print('relative residual norm: %e' % (diff.norm()/acq_data.norm()))

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
