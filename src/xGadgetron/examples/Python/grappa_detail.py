'''
Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
See function grappa_detail.m for an example showing more of the 
workings and functionality of the SIRF code.

Pre-requisites:
 1) This MATLAB code needs to be able to access a listening gadgetron.
    On the Virtual Machine, gadgetron is installed and the user just needs
    to type 'gadgetron' in a terminal window.
    On standalone systems, the user will need to have installed ISMRMRD
    and gadgetron code.

 2) An input data file from a GRAPPA MRI acquisition in the ISMRMRD format.
    Example GRAPPA datasets:
    a) 'meas_MID00108_FID57249_test_2D_2x.dat' is 
       available from https://www.ccppetmr.ac.uk/downloads
       This is in the manufacturer's raw data format and needs to be
       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
       This executable is installed on the Virtual Machine.

    b) A simulated ISMRMRD h5 file is available as default

Usage:
  grappa_basic.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of $SRC_PATH/SIRF
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import math
import matplotlib.pyplot as plt
import os
import sys

# import engine module
exec('from p' + args['--engine'] + ' import *')


def show(image_matrix, tile_shape, scale, titles):
    assert numpy.prod(tile_shape) >= image_matrix.shape[0],\
            "image tile rows x columns must equal the 3rd dim"\
            " extent of image_matrix"
    cols, rows = tile_shape
    vmin, vmax = scale
    fig = plt.figure()
    for z in range(image_matrix.shape[0]):
        ax = fig.add_subplot(cols, rows, z+1)
        ax.set_title(titles[z])
        ax.set_axis_off()
        imgplot = ax.imshow(image_matrix[z,:,:], vmin=vmin, vmax=vmax, cmap='gray')
    print('close figure 1 to continue')
    plt.show()


def main():
    
    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_filepath(data_path, args['--file'])
    
    # Initially we create a container that points to the h5 file. Data is
    # not read from file until the gadgetron is called using
    # the 'process' method.
    
    # Create an Acquisition Container of type pGadgetron.AcquisitionData
    print('---\n reading in file %s...' % input_file)
    input_Cont = AcquisitionData(input_file)
    
    
    # Pre-process this input data using three preparation gadgets
    # from gadgetron.
    # List gadgets to use (not all may be required for this test data).
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoAdjustROGadget', \
                                'RemoveROOversamplingGadget' ];
    
    # Call gadgetron by using the 'process' method. This runs the gadgets
    # specified in prep_gadgets, returning an instance
    # of an mGadgetron.AcquisitionsContainer
    preprocessed_AcCont = input_Cont.process(prep_gadgets);
    
    # Extract sorted k-space, permute dimensions and display
    ksp = preprocessed_AcCont.as_array(0);
    [ns,nc,nro] = preprocessed_AcCont.dimensions() ; # [nx ncoil ny]
    ksp = numpy.transpose(ksp,(1,2,0));
    show(abs(ksp[0,None,:,:]), tile_shape = (1,1), scale = (0, 0.7),\
            titles = ['Abs(Coil1)'])
            
    
    # Perform reconstruction of the preprocessed data.
    
    # 1) Create a recon object for the desired reconstruction.
    
    # In this demo, the recon object is created using the class
    # ImagesReconstructor(). A simpler class is available in the SIRF code
    # for a GRAPPA reconstruction:
    #   recon = GenericCartesianGRAPPAReconstruction()
    #
    #    To find what this does behind the scenes:
    #     type edit pGadgetron.GenericCartesianGRAPPAReconstruction
    #     and note the name assigned in the self function, here
    #       'SimpleGRAPPAReconstructionProcessor'.
    #     Then find the gadget chain defined by the class with the same
    #     name in the file xGadgetron/cGadgetron/chain_lib.h
    #
    
    recon_gadgets =  ['AcquisitionAccumulateTriggerGadget',
        'BucketToBufferGadget', 
        'GenericReconCartesianReferencePrepGadget', 
        'GRAPPA:GenericReconCartesianGrappaGadget', 
        'GenericReconFieldOfViewAdjustmentGadget', 
        'GenericReconImageArrayScalingGadget', 
        'ImageArraySplitGadget',
        ];
    
    recon = ImagesReconstructor(recon_gadgets) ;
    
    
    # 2) The GRAPPA gadget can compute G-factors in addition to
    # reconstructed images. We can set a gadget property as below if the gadget
    # has been identified with a label. In the above list of recon_gadgets,
    # the 4th is labelled 'GRAPPA' and we can use this label as below:
    recon.set_gadget_property('GRAPPA', 'send_out_gfactor', True)
    
    # If the chain had been set using
    # recon = GenericCartesianGRAPPAReconstruction(), an alternative method
    # would be available:
    #  recon.compute_gfactors(True);
    
    
    # 3) set the reconstruction input to be the data we just preprocessed.
    recon.set_input(preprocessed_AcCont);
    
    # 4) Run the reconstruction using 'process' to call gadgetron.
    print('---\n reconstructing...\n');
    recon.process();
    
    # Output
    
    # Reconstructed data sits in memory. We need to first get containers
    # for both the reconstructed images and g-factors, before extracting the
    # data as MATLAB arrays. Containers in effect point to the data.
    
    # Get images and gfactors as containers with type mGadgetron.ImagesContainer
    # (Note this syntax may change in the future with the addition of a
    #  method '.get_gfactor'.)
    images_imcont = recon.get_output('image');
    gfacts_imcont = recon.get_output('gfactor');
    
    # Return as MATLAB matrices the data pointed to by the containers.
    # Note the image data is complex.
    idata = images_imcont.as_array();
    maxv = numpy.amax(abs(idata))*0.6
    show(abs(idata[0,None,:,:]), tile_shape = (1,1), scale = (0, maxv),\
            titles = ['Abs(Image)'])
            
            
    gdata = gfacts_imcont.as_array();
    maxv = numpy.amax(abs(gdata))
    show(abs(gdata[0,None,:,:]), tile_shape = (1,1), scale = (0, maxv),\
            titles = ['G-factor map'])

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)



