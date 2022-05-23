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


data_path = args['<data_path>']


def main():
 
    # direct all engine's information and warnings printing to files
    msg_red = MessageRedirector('info.txt', 'warn.txt')

    cache_path = data_path;
    sens_filename = cache_path + "/sens_0.hv";
    print(sens_filename)
    tmpl_projdata_filename = data_path + "/tmpl_scanner.hs";
    print(tmpl_projdata_filename)

    acq_tmpl = AcquisitionData(tmpl_projdata_filename)
    print(acq_tmpl.dimensions())

    img_data = ImageData(sens_filename)
    print(img_data.dimensions())

    obj_fun = PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin()
    obj_fun.set_cache_path(path=cache_path, with_additive_corrections=True)

    print(obj_fun.get_cache_path())

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
