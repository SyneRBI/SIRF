'''sirf.Gadgetron Test set 4.
v{version}

2D Cartesian MR image reconstruction by direct creating
and running a chain of Gadgetron gadgets.

Usage:
  test4.py [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
'''

__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov"

import time
import sys

# import SIRF utilities
from sirf.Utilities import examples_data_path, existing_filepath, error
from sirf.Utilities import runner
# import MR engine types
from sirf.Gadgetron import AcquisitionData, Reconstructor

def test_main(rec=False, verb=False, throw=True):

    data_file = 'simulated_MR_2D_cartesian.h5'
    data_path = examples_data_path('MR')
    output_file = None
    type_to_save = 'all'
    show_plot = False
    algorithm = 'SimpleReconGadget'

    # locate the input data
    input_file = existing_filepath(data_path, data_file)
    acq_data = AcquisitionData(input_file)
    
    if algorithm == 'SimpleReconGadget':
        extra_gadgets = [algorithm]
    else:
        extra_gadgets = [algorithm, 'GenericReconFieldOfViewAdjustmentGadget']
    
    # create reconstruction object
    # Rather than using a predefined image reconstruction object, here a new 
    # image reconstruction object is created by concatinating multiple gadgets 
    # (for more information on Gadgetron and its gadgets please see: 
    # https://github.com/gadgetron/.).
    # Parameters for individual gadgets can be defined either during the 
    # creation of the reconstruction object:
    #   e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
    # or by giving a gadget a label (cf. label ex: for the last gadget)
    # and using set_gadget_property(label, propery, value).
    # The gadgets will be concatenated and will be executed as soon as 
    # process() is called.
    recon_gadgets = ['NoiseAdjustGadget',
        'AsymmetricEchoAdjustROGadget',
        'RemoveROOversamplingGadget',
        'AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)',
        'BucketToBufferGadget(split_slices=true, verbose=false)'] \
        + extra_gadgets + \
        ['ImageArraySplitGadget', 
        'ex:ExtractGadget'
        ]

    recon = Reconstructor(recon_gadgets)

    # ExtractGadget defines which type of image should be returned:
    # none      0
    # magnitude 1
    # real      2
    # imag      4
    # phase     8
    # in this example '5' returns both magnitude and imaginary part
##    recon.set_gadget_property('ex', 'extract_mask', 5)
    # === THE ABOVE IS OBSOLETE, NOW SHOULD USE ===>
    if type_to_save=='mag' or type_to_save=='all':
        recon.set_gadget_property('ex', 'extract_magnitude', True)
    if type_to_save=='imag' or type_to_save=='all':
        recon.set_gadget_property('ex', 'extract_imag', True)
    
    # provide raw k-space data as input
    recon.set_input(acq_data)

    # optionally set Gadgetron server host and port
    recon.set_host('localhost')
    # On VM you can try a port other than the default 9002, e.g. 9003, by taking
    # the following steps:
    # 1) in ~/devel/install/share/gadgetron/config/gadgetron.xml replace
    #    <port>9002</port> with <port>9003</port>
    # 2) go to Settings->Network->Advanced->Port Forwarding and add new rule
    #    (click on green + in the upper right corner) with Host and Guest ports
    #    set to 9003
    # 3) uncomment the next line
    #recon.set_port('9003')
    # Note: each gadget chain can run on a different VM - to try, start two VMs
    # and do the above steps 1 and 2 on one of them, then add
    # recon.set_port('9003') before recon.process in grappa_detail.py
    # (where preprocessing will still run on default port 9002). 

    # perform reconstruction
    recon.process()
    
    # retrieve reconstructed image data
    image_data = recon.get_output()

    # show reconstructed image data
    if show_plot:
        for im in range(image_data.number()):
            image = image_data.image(im)
            # image types   series
            # magnitude 1       0
            # phase     2    3000
            # real      3    1000
            # imag      4    2000
            im_type = image.image_type()
            im_series = image.image_series_index()
            print('image: %d, type: %d, series: %d' % (im, im_type, im_series))
        image_data.show(title = 'Images magnitude and imaginary part')

    if output_file is not None:
        filename = output_file
        i = filename.find('.')
        if i < 0:
            ext = 'h5'
        else:
            ext = filename[i + 1:]
            filename = filename[:i]
        print('writing to %s' % (filename + '.' + ext))
        image_data.write(filename, ext=ext)

    return 0, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
