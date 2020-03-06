'''OSEM reconstruction demo.
We actually use the OSMAPOSL reconstructor in this demo. This reconstructor
implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
from Green et al for Maximum a Posteriori (MAP) maximisation. Here we use it
for Maximum Likelihood (ML) in which case it is equivalent to OSEM.

Usage:
  osem_reconstruction [--help | options]

Options:
  -S <file>, --sino=<file>     sinogram [default: data/examples/PET/mMR/my_forward_projection.hs]
  -R <file>, --rand=<file>     randoms
  -a <attn>, --attn=<attn>     attenuation image file
  -n <norm>, --norm=<norm>     ECAT8 bin normalization file
  -s <subs>, --subs=<subs>     number of subsets [default: 12]
  -i <siter>, --subiter=<siter>  number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -o <outp>, --outp=<outp>     output file prefix [default: recon]
  -d <nxny>, --nxny=<nxny>     image x and y dimensions as string '(nx,ny)'
                               (no space after comma) [default: (127,127)]
  --visualisations             show visualisations
  --nifti                      save output as nifti
  --gpu                        use gpu
'''

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

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from ast import literal_eval


def file_exists(filename):
    """Check if file exists, optionally throw error if not"""
    return os.path.isfile(filename)


def check_file_exists(filename):
    """Check file exists, else throw error"""
    if not file_exists:
        raise error('File not found: %s' % filename)


# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

# Sinogram. if sino not found, get the one in the example data
sino_file = args['--sino']
if not file_exists(sino_file):
    sino_file = examples_data_path('PET') + '/mMR/' + sino_file
    check_file_exists(sino_file)

# Randoms
if args['--rand']:
    rand_file = args['--rand']
    check_file_exists(rand_file)

# Attenuation
if args['--attn']:
    attn_file = args['--attn']
    check_file_exists(attn_file)

# Norm
if args['--norm']:
    norm_file = args['--norm']
    check_file_exists(norm_file)

# Number of voxels
nxny = literal_eval(args['--nxny'])

# Output file
outp_file = args['--outp']

if args['--visualisations']:
    visualisations = True
else:
    visualisations = False

if args['--nifti']:
    nifti = True
    import sirf.Reg
else:
    nifti = False

if args['--gpu']:
    use_gpu = True
else:
    use_gpu = False

# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])


def main():
 
    # direct all engine's information and warnings printing to files
    msg_red = MessageRedirector('info.txt', 'warn.txt')

    if not use_gpu:
        print("Using CPU projector...")
        # select acquisition model that implements the geometric
        # forward projection by a ray tracing matrix multiplication
        acq_model = AcquisitionModelUsingRayTracingMatrix()
    else:
        print("Using GPU projector...")
        # If using GPU, use the niftypet projector
        acq_model = AcquisitionModelUsingNiftyPET()

    # PET acquisition data to be read from this file
    # (TODO: a link to raw data formats document to be given here)
    print('raw data: %s' % sino_file)
    acq_data = AcquisitionData(sino_file)

    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object ad) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0, nxny)
    # If using GPU, need to make sure that image is right size.
    if use_gpu:
        image.initialise(dim=(127,320,320), vsize=(2.03125,2.08626,2.08626))
        image.fill(1.0)

    # If norm is present
    if norm_file:
        # create acquisition sensitivity model from ECAT8 normalisation data
        asm_norm = AcquisitionSensitivityModel(norm_file)
    
    # If attenuation is present
    if attn_file:
        attn_image = ImageData(attn_file)
        # If gpu, make sure that the attn. image is same dimensions as image to recon
        if use_gpu:
            resampler = sirf.Reg.NiftyResample()
            resampler.set_reference_image(image)
            resampler.set_floating_image(attn_image)
            resampler.set_interpolation_type_to_linear()
            resampler.set_padding_value(0.0)
            attn_image = resampler.forward(attn_image)
        asm_attn = AcquisitionSensitivityModel(attn_image, acq_model)
        # temporary fix pending attenuation offset fix in STIR:
        # converting attenuation into 'bin efficiency'
        asm_attn.set_up(acq_data)
        bin_eff = AcquisitionData(acq_data)
        bin_eff.fill(1.0)
        print('applying attenuation (please wait, may take a while)...')
        asm_attn.unnormalise(bin_eff)
        asm_attn = AcquisitionSensitivityModel(bin_eff)

    # Get ASM dependent on attn and/or norm
    if norm_file and attn_file:
        print("AcquisitionSensitivityModel contains norm and attenuation...")
        asm = AcquisitionSensitivityModel(asm_norm, asm_attn)
    elif norm_file:
        print("AcquisitionSensitivityModel contains norm...")
        asm = asm_norm
    elif attn_file:
        print("AcquisitionSensitivityModel contains attenuation...")
        asm = asm_attn
    if asm:
        print("Setting AcquisitionSensitivityModel...")
        acq_model.set_acquisition_sensitivity(asm)

    # If randoms are present
    if rand_file:
        print("Adding randoms...")
        randoms = AcquisitionData(rand_file)
        acq_model.set_background_term(randoms)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm (since we are not using a penalty, or prior, in
    # this example, we actually run OSEM);
    # this algorithm does not converge to the maximum of the objective function
    # but is used in practice to speed-up calculations
    recon = OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)
    recon.set_input(acq_data)

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
        sirf.Reg.NiftiImageData(out).write(outp_file)

    if visualisations:
        # show reconstructed image
        image_array = out.as_array()
        show_2D_array('Reconstructed image', image_array[z,:,:])
        pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
