'''Listmode reconstruction demo

Usage:
  listmode_reconstruction [--help | options]

Options:
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET/mMR
                               subfolder of SIRF root folder
  -l <list>, --list=<list>     listmode file [default: list.l.hdr]
  -g <sino>, --sino=<sino>     output file prefix [default: sinograms]
  -a <attn>, --attn=<attn>     attenuation image file file [default: mu_map.hv]
  -n <norm>, --norm=<norm>     ECAT8 bin normalization file [default: norm.n.hdr]
  -I <int>, --interval=<int>   scanning time interval to convert as string '(a,b)'
                               (no space after comma) [default: (0,50)]
  -d <nxny>, --nxny=<nxny>     image x and y dimensions as string '(nx,ny)'
                               (no space after comma) [default: (127,127)]
  -S <subs>, --subs=<subs>     number of subsets [default: 7]
  -i <siter>, --subiter=<siter>  number of sub-iterations [default: 2]
  -o <outp>, --outp=<outp>     output file prefix [default: recon]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
  -C <cnts>, --counts=<cnts>   account for delay between injection and acquisition start by shifting interval to start when counts exceed given threshold.
  --visualisations             show visualisations
  --nifti                      save output as nifti
  --gpu                        use gpu
  --non-interactive            do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2018 - 2021, 2024 University College London.
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

__version__ = '1.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from ast import literal_eval
import os

import PET_plot_functions

from sirf.Utilities import error, examples_data_path, existing_filepath
from sirf.Utilities import show_2D_array

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)


# process command-line options
data_path = args['--path']
if data_path is None:
    # default to data/examples/PET/mMR
    # Note: seem to need / even on Windows
    #data_path = os.path.join(examples_data_path('PET'), 'mMR')
    data_path = examples_data_path('PET') + '/mMR'
print('Finding files in %s' % data_path)

list_file = args['--list']
sino_file = args['--sino']
tmpl_file = args['--tmpl']
norm_file = args['--norm']
attn_file = args['--attn']
outp_file = args['--outp']
# Check file exists (e.g., absolute path). Else prepend data_path
if not os.path.isfile(list_file):
    list_file = existing_filepath(data_path, list_file)
if not os.path.isfile(tmpl_file):
    tmpl_file = existing_filepath(data_path, tmpl_file)
if not os.path.isfile(norm_file):
    norm_file = existing_filepath(data_path, norm_file)
if not os.path.isfile(attn_file):
    attn_file = existing_filepath(data_path, attn_file)
nxny = literal_eval(args['--nxny'])
input_interval = literal_eval(args['--interval'])
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
storage = args['--storage']
count_threshold = args['--counts']

if args['--visualisations']:
    visualisations = True
else:
    visualisations = False
if args['--non-interactive']:
    visualisations = False

if args['--gpu']:
    use_gpu = True
#    import sirf.Reg
else:
    use_gpu = False


def main():

    # engine's messages go to files, except error messages, which go to stdout
    _ = pet.MessageRedirector('info.txt', 'warn.txt')

    # select acquisition data storage scheme
    pet.AcquisitionData.set_storage_scheme(storage)

    listmode_data = pet.ListmodeData(list_file)

    # First step is to create AcquisitionData ("sinograms") from the
    # listmode file.
    # See the listmode_to_sinograms demo for some more information on this step.

    # create listmode-to-sinograms converter object
    # See also the listmode_to_sinograms demo
    lm2sino = pet.ListmodeToSinograms()

    # set input, output and template files
    lm2sino.set_input(listmode_data)
    lm2sino.set_output_prefix(sino_file)
    # need to be at maximum resolution to work in listmode reconstruction
    acq_data_template = listmode_data.acquisition_data_template()
    #print(acq_data_template.get_info())
    lm2sino.set_template(acq_data_template)

    if count_threshold is None:
        interval = input_interval
    else:
        time_shift = lm2sino.get_time_at_which_num_prompts_exceeds_threshold(count_threshold)
        if time_shift < 0:
            print("No time found at which count rate exceeds " + str(time_shift) + ", not modifying interval")
        interval = (input_interval[0]+time_shift, input_interval[1]+time_shift)
        print("Time at which count rate exceeds " + str(count_threshold) + " = " + str(time_shift) + " s.")
        print("Input intervals: " + str(input_interval[0]) + ", " + str(input_interval[1]))
        print("Modified intervals: " + str(interval[0]) + ", " + str(interval[1]))

    # set interval
    lm2sino.set_time_interval(interval[0], interval[1])

    # set up the converter
    lm2sino.set_up()

    # convert (need it for the scatter estimate)
    lm2sino.process()

    acq_data = lm2sino.get_output()

    # Get the randoms
    randoms = lm2sino.estimate_randoms()

    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object acq_data) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0, nxny)

    # read attenuation image
    attn_image = pet.ImageData(attn_file)
    if visualisations:
        z = attn_image.shape[0]//2
        attn_image_as_array = attn_image.as_array()
        show_2D_array('Attenuation image', attn_image_as_array[z,:,:])

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()
    acq_model.set_num_tangential_LORs(10)

    # create acquisition sensitivity model from normalisation data
    asm_norm = pet.AcquisitionSensitivityModel(norm_file)

    # create attenuation factors
    asm_attn = pet.AcquisitionSensitivityModel(attn_image, acq_model)
    # converting attenuation image into attenuation factors (one for every bin)
    asm_attn.set_up(acq_data)
    ac_factors = acq_data.get_uniform_copy(value=1)
    print('applying attenuation (please wait, may take a while)...')
    asm_attn.unnormalise(ac_factors)
    asm_attn = pet.AcquisitionSensitivityModel(ac_factors)

    # scatter estimation
    print('estimating scatter (this will take a while!)')
    scatter_estimator = pet.ScatterEstimator()

    scatter_estimator.set_input(acq_data)
    scatter_estimator.set_attenuation_image(attn_image)
    scatter_estimator.set_randoms(randoms)
    scatter_estimator.set_asm(asm_norm)
    # invert attenuation factors to get the correction factors,
    # as this is unfortunately what a ScatterEstimator needs
    acf_factors=acq_data.get_uniform_copy()
    acf_factors.fill(1/ac_factors.as_array())
    scatter_estimator.set_attenuation_correction_factors(acf_factors)
    scatter_estimator.set_output_prefix(sino_file + '_scatter')
    scatter_estimator.set_num_iterations(3)
    scatter_estimator.set_up()
    scatter_estimator.process()
    scatter_estimate = scatter_estimator.get_output()
    if visualisations:
        scatter_estimate_as_array = scatter_estimate.as_array()
        show_2D_array('Scatter estimate (first sinogram)', scatter_estimate_as_array[0, 0, :, :])
        PET_plot_functions.plot_sinogram_profile(acq_data, randoms=randoms, scatter=scatter_estimate)

    # chain attenuation and ECAT8 normalisation
    asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    asm.set_up(acq_data)

    acq_model.set_acquisition_sensitivity(asm)
    acq_model.set_background_term(randoms + scatter_estimate)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = pet.PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin();
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_acquisition_data(listmode_data)

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm (since we are not using a penalty, or prior, in
    # this example, we actually run OSEM);
    # this algorithm does not converge to the maximum of the objective function
    # but is used in practice to speed-up calculations
    # See the reconstruction demos for more complicated examples
    recon = pet.OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)

    # set up the reconstructor based on a sample image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    recon.set_up(image)

    # set the initial image estimate
    recon.set_current_estimate(image)

    # reconstruct
    print('reconstructing, please wait...')
    recon.process()
    out = recon.get_output()
    if not args['--nifti']:
        out.write(outp_file)
    else:
        import sirf.Reg as reg
        reg.NiftiImageData(out).write(outp_file)

    if visualisations:
        # show reconstructed image
        z = out.shape[0]//2
        image_array = out.as_array()
        show_2D_array('Reconstructed image', image_array[z,:,:])
#        pylab.show()


try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    print('%s' % err.value)
