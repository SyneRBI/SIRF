'''OSEM reconstruction demo.
We actually use the OSMAPOSL reconstructor in this demo. This reconstructor
implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
from Green et al for Maximum a Posteriori (MAP) maximisation. Here we use it
for Maximum Likelihood (ML) in which case it is equivalent to OSEM.

Usage:
  osem_reconstruction [--help | options]

Options:
  -S <file>, --sino=<file>     sinogram (default: data/examples/PET/my_forward_projection.hs)
  -R <file>, --rand=<file>     randoms
  -a <attn>, --attn_im=<attn>  attenuation image file
  -A <attn>, --attn_sn=<attn>  attenuation sinogram
  -n <norm>, --norm_e8=<norm>  ECAT8 bin normalization file
  -N <norm>, --norm_sn=<norm>  Normalisation sinogram
  -s <subs>, --subs=<subs>     number of subsets [default: 12]
  -i <siter>, --subiter=<siter>  number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -o <outp>, --outp=<outp>     output file prefix [default: recon]
  -d <nxny>, --nxny=<nxny>     image x and y dimensions as string '(nx,ny)'
                               (no space after comma) [default: (127,127)]
  --visualisations             show visualisations
  --nifti                      save output as nifti
  --gpu                        use gpu
  --non-interactive            do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2020 University College London.
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
if not sino_file:
    print("Sinogram not given, using data/examples/PET/my_forward_projection.hs")
    sino_file = examples_data_path('PET') + "/my_forward_projection.hs"
if not file_exists(sino_file):
    raise error("Sinogram not found: " + sino_file)

# Randoms
rand_file = None
if args['--rand']:
    rand_file = args['--rand']
    check_file_exists(rand_file)

# Attenuation - image
attn_im_file = None
if args['--attn_im']:
    attn_im_file = args['--attn_im']
    check_file_exists(attn_im_file)

# Attenuation - sinogram
attn_sn_file = None
if args['--attn_sn']:
    attn_sn_file = args['--attn_sn']
    check_file_exists(attn_sn_file)

# Norm - ECAT8
norm_e8_file = None
if args['--norm_e8']:
    norm_e8_file = args['--norm_e8']
    check_file_exists(norm_e8_file)

# Norm - sinogram
norm_sn_file = None
if args['--norm_sn']:
    norm_sn_file = args['--norm_sn']
    check_file_exists(norm_sn_file)

# Number of voxels
nxny = literal_eval(args['--nxny'])

# Output file
outp_file = args['--outp']

if args['--visualisations']:
    visualisations = True
else:
    visualisations = False
if args['--non-interactive']:
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

    if not use_gpu:
        print("Using CPU projector...")
        # select acquisition model that implements the geometric
        # forward projection by a ray tracing matrix multiplication
        acq_model = AcquisitionModelUsingRayTracingMatrix()
    else:
        print("Using GPU projector...")
        # If using GPU, use the niftypet projector
        acq_model = AcquisitionModelUsingNiftyPET()
        # Truncate to FOV
        acq_model.set_use_truncation(True)

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
    asm_norm = None
    if norm_e8_file and norm_sn_file:
        raise error("For normalisation, only give ECAT8 or sinogram.")
    if norm_e8_file:
        # create acquisition sensitivity model from ECAT8 normalisation data
        asm_norm = AcquisitionSensitivityModel(norm_e8_file)
    if norm_sn_file:
        norm_sino = AcquisitionData(norm_sn_file)
        asm_norm = AcquisitionSensitivityModel(norm_sino)
    
    # If attenuation is present
    asm_attn = None
    if attn_im_file and attn_sn_file:
        raise error("For attenuation, only give image or sinogram.")
    if attn_im_file:
        attn_image = ImageData(attn_im_file)
        # If gpu, make sure that the attn. image is same dimensions as image to recon
        if use_gpu:
            resampler = sirf.Reg.NiftyResampler()
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
    if attn_sn_file:
        attn_sino = AcquisitionData(attn_sn_file)
        asm_attn = AcquisitionSensitivityModel(attn_sino)

    # Get ASM dependent on attn and/or norm
    asm = None
    if asm_norm and asm_attn:
        print("AcquisitionSensitivityModel contains norm and attenuation...")
        asm = AcquisitionSensitivityModel(asm_norm, asm_attn)
    elif asm_norm:
        print("AcquisitionSensitivityModel contains norm...")
        asm = asm_norm
    elif asm_attn:
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
        z = image_array.shape[0]//3
        show_2D_array('Reconstructed image', image_array[z,:,:])
        pylab.show()


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
