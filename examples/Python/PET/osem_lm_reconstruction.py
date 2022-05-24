'''OSEM reconstruction from listmode demo.


Usage:
  osem_lm_reconstruction [--help | options] <data_path>
  
Options:
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
'''


__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)


# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


data_path = args['<data_path>'] + '/'


def main():
 
    # direct all engine's information and warnings printing to files
    msg_red = MessageRedirector('info.txt', 'warn.txt')

    cache_path = data_path;
    sens_filename = cache_path + "sens_0.hv";
    #print(sens_filename)
    tmpl_projdata_filename = data_path + "tmpl_scanner.hs";
    #print(tmpl_projdata_filename)

    acq_tmpl = AcquisitionData(tmpl_projdata_filename)
    #print(acq_tmpl.dimensions())

    img_data = ImageData(sens_filename)
    #print(img_data.dimensions())

    #obj_fun = PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin()
    obj_fun = make_Poisson_loglikelihood( \
        likelihood_type='LinearModelForMeanAndListModeDataWithProjMatrixByBin')
    obj_fun.set_cache_path(cache_path, True)
    obj_fun.set_skip_lm_input_file(True);
    obj_fun.set_skip_balanced_subsets(True);
    obj_fun.set_acquisition_data(acq_tmpl);
    obj_fun.set_skip_balanced_subsets(True);
    obj_fun.set_max_ring_difference(60);
    obj_fun.set_cache_max_size(1500000000);

    #print(obj_fun.get_cache_path())
    #print(obj_fun.get_cache_max_size())
    #print(obj_fun.get_subsensitivity_filenames())

    num_subiterations = 11
    recon = OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subiterations(num_subiterations)
    recon.set_output_filename_prefix(data_path + "/my_output");
    recon.set_save_interval(1);

    print('setting up, please wait...')
    recon.set_up(img_data)

    print('reconstructing, please wait...')
    recon.reconstruct(img_data)

    img_out = recon.get_output()
    print(img_out.dimensions())

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
