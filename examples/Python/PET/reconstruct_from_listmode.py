'''Listmode-to-sinograms conversion and reconstruction demo.

Usage:
  reconstruct_from_listmode [--help | options]

Options:
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -l <list>, --list=<list>     listmode file [default: list.l.hdr.STIR]
  -g <sino>, --sino=<sino>     output file prefix [default: sinograms]
  -t <tmpl>, --tmpl=<tmpl>     raw data template [default: template_span11.hs]
  -a <attn>, --attn=<attn>     attenuation image file file [default: mu_map.hv]
  -n <norm>, --norm=<norm>     ECAT8 bin normalization file [default: norm.n.hdr.STIR]
  -i <int>, --interval=<int>   scanning time interval to convert as string '(a,b)'
                               (no space after comma) [default: (0,100)]
  -d <nxny>, --nxny=<nxny>     image x and y dimensions as string '(nx,ny)'
                               (no space after comma) [default: (127,127)]
  -S <subs>, --subs=<subs>     number of subsets [default: 2]
  -I <iter>, --iter=<iter>     number of iterations [default: 2]
  -o <outp>, --outp=<outp>     output file prefix [default: recon]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
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

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
list_file = args['--list']
sino_file = args['--sino']
tmpl_file = args['--tmpl']
norm_file = args['--norm']
attn_file = args['--attn']
outp_file = args['--outp']
list_file = existing_filepath(data_path, list_file)
tmpl_file = existing_filepath(data_path, tmpl_file)
norm_file = existing_filepath(data_path, norm_file)
attn_file = existing_filepath(data_path, attn_file)
nxny = literal_eval(args['--nxny'])
interval = literal_eval(args['--interval'])
num_subsets = int(args['--subs'])
num_iterations = int(args['--iter'])
storage = args['--storage']

def main():

    # output goes to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # create listmode-to-sinograms converter object
    lm2sino = ListmodeToSinograms()

    # set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix(sino_file)
    lm2sino.set_template(tmpl_file)

    # set interval
    lm2sino.set_time_interval(interval[0], interval[1])

    # switch flags on/off
    lm2sino.flag_on('store_prompts')
    lm2sino.flag_off('interactive')

    # set up the converter
    lm2sino.set_up()

    # convert
    lm2sino.process()

    # get access to the sinograms
    acq_data = lm2sino.get_output()
    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions: %dx%dx%d' % acq_dim)
    z = acq_dim[0]//2
    show_2D_array('Acquisition data', acq_array[z,:,:])

    # read attenuation image
    attn_image = ImageData(attn_file)
    attn_image_as_array = attn_image.as_array()
    z = attn_image_as_array.shape[0]//2
    show_2D_array('Attenuation image', attn_image_as_array[z,:,:])

    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object ad) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0, nxny)

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix()
    acq_model.set_up(acq_data, attn_image)

    # create acquisition sensitivity model from ECAT8 normalization data
    asm_norm = AcquisitionSensitivityModel(norm_file)

    asm_attn = AcquisitionSensitivityModel(attn_image, acq_model)
    asm_attn.set_up(acq_data)
    bin_eff = AcquisitionData(acq_data)
    bin_eff.fill(1.0)
    print('applying attenuation (please wait, may take a while)...')
    asm_attn.unnormalise(bin_eff)
    asm_beff = AcquisitionSensitivityModel(bin_eff)

    asm = AcquisitionSensitivityModel(asm_norm, asm_beff)
    acq_model.set_acquisition_sensitivity(asm)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_max_segment_num_to_process(1)

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm (since we are not using a penalty, or prior, in
    # this example, we actually run OSEM);
    # this algorithm does not converge to the maximum of the objective function
    # but is used in practice to speed-up calculations
    recon = OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_iterations)

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

    # show reconstructed image
    image_array = recon.get_current_estimate().as_array()
    show_2D_array('Reconstructed image', image_array[20,:,:])
    pylab.show()

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
