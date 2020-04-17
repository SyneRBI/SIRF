'''Regsitration of SIRF images.

Usage:
  registration [--algo=<algo>] [--ref=<file>,<eng>] [--flo=<file>,<eng>]... [--warped_prefix=<file>]
  [--disp_fwd_prefix=<file>] [--def_fwd_prefix=<file>] [--disp_inv_prefix=<file>] [--def_inv_prefix=<file>]
  [--TM_fwd_prefix=<file>] [--TM_inv_prefix=<file>] [--rmask=<file>,<eng>] [--fmask=<file>,<eng>] [--print]
  [--par_file=<file>] [--par=<arg>]... [--working_folder=<file>] [--overwrite=<bool>] [--delete=<bool>]
  [--help | options]

Options:
  --ref <file>,<eng>           reference image (default: test.nii.gz,Reg)
  --flo <file>,<eng>           floating image (default: test2.nii.gz,Reg)
  --algo <algo>                registration algorithm (aladin,f3d,spm) [default: aladin]

  --warped_prefix <file>       warped image filename prefix
  --disp_fwd_prefix <file>     forward displacement field image
  --def_fwd_prefix <file>      forward deformation field image
  --disp_inv_prefix <file>     inverse displacement field image
  --def_inv_prefix <file>      inverse deformation field image

  --TM_fwd_prefix <file>       forward transformation matrix (rigid/affine only)
  --TM_inv_prefix <file>       inverse transformation matrix (rigid/affine only)

  --rmask <file,eng>           mask of reference image (NiftyReg only)
  --fmask <file,eng>           mask of floating image (NiftyReg only)
  --print                      print wrapped methods (NiftyReg only)
  --par_file <file>            parameter file (NiftyReg only)
  --par <string>               set wrapped parameter (NiftyReg only). Examples: '--par SetPerformAffine,0', '--par SetInterpolationToCubic', '--par SetFloatingThresholdUp,1,2')

  --working_folder <fname>     folder in which to save temporary files (SPM only) (default: cwd/spm_working_folder)
  --overwrite <bool>           should I overwrite files if already present? (SPM only) (default: true)
  --delete <bool>              should I delete temporary files? (SPM only) (default: True)
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2020 University College London.
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

import os
import sirf.Reg
from sirf.Utilities import *

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

def get_image(arg):
    """Get an image filename and its engine"""
    arg = arg.split(',')
    if len(arg) != 2:
        raise error('image argument expects filename and engine')
    filename = arg[0]
    engine = arg[1]
    if engine == 'Reg':
        global sirf
        im = sirf.Reg.ImageData(filename)
    elif engine == 'STIR':
        import sirf.STIR
        im = sirf.STIR.ImageData(filename)
    elif engine == 'Gadgetron':
        import sirf.Gadgetron
        im = sirf.Gadgetron.ImageData(filename)
    else:
        raise error('unknown engine: ' + engine)
    return im

def get_algorithm(algo):
    """Get algorithm based on string"""
    if algo == "aladin":
        reg = sirf.Reg.NiftyAladinSym()
    elif algo == "f3d":
        reg = sirf.Reg.NiftyF3dSym()
    elif algo == "spm":
        reg = sirf.Reg.SPMRegistration()
    else:
        raise error('unknown algorithm')
    return reg


def main():

    # files
    ref_args = args['--ref']
    flo_args = args['--flo']
    # if using the default for any, need to get the examples folder
    if ref_args is None or flo_args is None:
        SIRF_PATH = os.environ.get('SIRF_PATH')
        if SIRF_PATH is not None:
            examples_path = SIRF_PATH + '/data/examples/Registration'
        else:
            errorMsg = 'You need to set the SIRF_PATH environment variable to allow finding the raw data.'
            raise error(errorMsg)
    # Ref
    if ref_args is None:
        ref_args = examples_path + '/test.nii.gz,Reg'
    ref = get_image(ref_args)
    # Flo
    if len(flo_args) == 0:
        flo_args = [examples_path + '/test2.nii.gz,Reg']
    flos = []
    for flo_arg in flo_args:
            flos.append(get_image(flo_arg))

    # Algorithm
    algo = args['--algo']
    reg = get_algorithm(algo)

    if args['--print']:
        if algo == 'spm':
            raise error('--print only available for NiftyReg')
        reg.print_all_wrapped_methods()
        exit(0)

    # Set images
    reg.set_reference_image(ref)
    for flo in flos:
        reg.add_floating_image(flo)

    # rmask
    if args['--rmask']:
        if algo == 'spm':
            raise error('--rmask only available for NiftyReg')
        reg.set_reference_mask(get_image(args['--rmask']))
    # fmask
    if args['--fmask']:
        if algo == 'spm':
            raise error('--fmask only available for NiftyReg')
        reg.set_floating_mask(get_image(args['--fmask']))
    # par_file
    if args['--par_file']:
        if algo == 'spm':
            raise error('--par_file only available for NiftyReg')
        reg.set_parameter_file(args['--par_file'])
    # pars
    if args['--par_file']:
        if algo == 'spm':
            raise error('--par_file only available for NiftyReg')
        reg.set_parameter_file(args['--par_file'])

    pars = args['--par']
    if len(pars) > 0:
        if algo == 'spm':
            raise error('--par only available for NiftyReg')
        for par in pars:
            par_split = par.split(',')
            if len(par_split) == 1:
                reg.set_parameter(par_split[0])
            elif len(par_split) == 2:
                reg.set_parameter(par_split[0], par_split[1])
            elif len(par_split) == 3:
                reg.set_parameter(par_split[0], par_split[1], par_split[2])
            else:
                raise error('Max number of NiftyReg args is 2.')

    # working folder
    working_folder = args['--working_folder']
    if algo == 'spm' and working_folder is None:
        working_folder = os.getcwd() + '/spm_working_folder'
    if working_folder is not None:
        if algo != 'spm':
            raise error('--working_folder only available for spm')
        reg.set_working_folder(working_folder)

    # overwrite
    overwrite = args['--overwrite']
    if algo == 'spm' and overwrite is None:
        overwrite = True
    if overwrite:
        if algo != 'spm':
            raise error('--overwrite only available for spm')
        reg.set_working_folder_file_overwrite(overwrite)

    # delete temp files
    delete_temp_files = args['--delete']
    if algo == 'spm' and delete_temp_files is None:
        delete_temp_files = True
    if delete_temp_files:
        if algo != 'spm':
            raise error('--delete only available for spm')
        reg.set_delete_temp_files(delete_temp_files)

    # Register
    reg.process()

    # Save results
    warped_prefix = args['--warped_prefix']
    disp_fwd_prefix = args['--disp_fwd_prefix']
    disp_inv_prefix = args['--disp_inv_prefix']
    def_fwd_prefix = args['--def_fwd_prefix']
    def_inv_prefix = args['--def_inv_prefix']
    TM_fwd_prefix = args['--TM_fwd_prefix']
    TM_inv_prefix = args['--TM_inv_prefix']
    for i in range(len(flos)):

        if warped_prefix is not None:
            reg.get_output(i).write(warped_prefix + str(i))
        if disp_fwd_prefix is not None:
            reg.get_displacement_field_forward(i).write(disp_fwd_prefix + str(i))
        if disp_inv_prefix is not None:
            reg.get_displacement_field_inverse(i).write(disp_inv_prefix + str(i))
        if def_fwd_prefix is not None:
            reg.get_deformation_field_forward(i).write(def_fwd_prefix + str(i))
        if def_inv_prefix is not None:
            reg.get_deformation_field_inverse(i).write(def_inv_prefix + str(i))
        if TM_fwd_prefix is not None:
            reg.get_transformation_matrix_forward(i).write(TM_fwd_prefix + str(i))
        if TM_inv_prefix is not None:
            reg.get_transformation_matrix_inverse(i).write(TM_inv_prefix + str(i))

main()
