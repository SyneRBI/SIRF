function grappa_detail
% GRAPPA_DETAIL Annotated demo for reconstruction of GRAPPA acquired data.
% 
% Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
% See function grappa_basic.m for a simpler example.
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file from a GRAPPA MRI acquisition in the ISMRMRD format.
%    Example GRAPPA datasets:
%    a) 'meas_MID00108_FID57249_test_2D_2x.dat' is 
%       available from https://www.ccppetmr.ac.uk/downloads
%       This is in the manufacturer's raw data format and needs to be
%       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
%       This executable is installed on the Virtual Machine.
%
%    b) An ISMRMRD h5 file can be simulated using the MATLAB function
%    gen_us_data.m . This will reconstruct faster than the real data in a).
%
% Usage:
%  grappa_detail
%
%
% Adapted by David Atkinson (D.Atkinson@ucl.ac.uk) from original 
% code by Evgueni Ovtchinnikov
%
% See also GRAPPA_BASIC GEN_US_DATA


% load the SIRF mutilities library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

% This demo uses gadgetron as the MR reconstruction engine. 
% Load the SIRF library and import the SIRF class.
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

import mGadgetron.*


%try  Commented out the try-catch block so that we get more informative
%     errors from Matlab.

    % Get the filename of the input ISMRMRD h5 file
    disp('Select ISMRMRD H5 file')
    [fn,pn] = uigetfile('*','Select ISMRMRD H5 file') ;
    filein = fullfile(pn,fn) ;
    
    % Initially we create a container that points to the h5 file. Data will
    % not be read from file until the gadgetron is called using 
    % the 'process' method.
    
    % Create an Acquisition Container. Here because of the previous 
    % 'import mGadgetron.*', this will be of type mGadgetron.AcquisitionData
    input_data = AcquisitionData(filein);
    
    % Pre-process this input data using three preparation gadgets 
    % from gadgetron.
    % List gadgets to use (not all are required for this test data).
    prep_gadgets = [{'NoiseAdjustGadget'}, ...
                    {'AsymmetricEchoAdjustROGadget'} ...
                    {'RemoveROOversamplingGadget'} ];
    
    % Call gadgetron by using the 'process' method. This runs the gadgets 
    % specified in prep_gadgets, returning an instance 
    % of an mGadgetron.AcquisitionsContainer
    preprocessed_AcCont = input_data.process(prep_gadgets);
    
    
    % Perform reconstruction of the preprocessed data.
    
    % The recon object can be created using the class 
    % ImagesReconstructor(), or for this example, 
    % GenericCartesianGRAPPAReconstruction(). Currently, not all methods 
    % in the 2nd class are defined in ImagesReconstructor() and so this 
    % demo is limited to using GenericCartesianGRAPPAReconstruction(). 
    
    % 1) Create a recon object for the desired reconstruction. 
    %    
    recon = GenericCartesianGRAPPAReconstruction();
    
    %
    %    This creates the recon object with type 
    %    mGadgetron.GenericCartesianGRAPPAReconstruction
    %    set up for a GRAPPA reconstruction.
    %
    %    To find what this does behind the scenes:
    %     type edit mGadgetron.GenericCartesianGRAPPAReconstruction
    %     and note the name assigned in the self function, here 
    %       'SimpleGRAPPAReconstructionProcessor'.
    %     Then find the gadget chain defined by the class with the same 
    %     name in the file xGadgetron/cGadgetron/chain_lib.h
    %
    %     In this case, it is a chain of 7 gadgets. In the pre-processing
    %     step above we explicitly chained gadgets and in principle we
    %     could do that here using the ImagesReconstructor():
%     
%     recon_gadgets =  [...
%         {'AcquisitionAccumulateTriggerGadget'}, ...
%         {'BucketToBufferGadget'}, ...
%         {'GenericReconCartesianReferencePrepGadget'}, ...
%         {'GRAPPA:GenericReconCartesianGrappaGadget'}, ...
%         {'GenericReconFieldOfViewAdjustmentGadget'}, ...
%         {'GenericReconImageArrayScalingGadget'}, ...
%         {'ImageArraySplitGadget'} ...
%         ];
%     
%     recon = ImagesReconstructor(recon_gadgets) ;
    %
    % However, due to limitations on outputing the gfactors, this demo
    % does not use the ImagesReconstructor() ;
    
    % 2) The GRAPPA gadget can compute G-factors in addition to 
    % reconstructed images. To set this property to 'true', if a gadget 
    % has been labelled, as in the 4th gadget above, we can 
    % type; 
    %   recon.set_gadget_property('GRAPPA', 'send_out_gfactor', false)
    % or if the chain was set using 
    % recon = GenericCartesianGRAPPAReconstruction() we can do:
       
    recon.compute_gfactors(true);
    %   NOTE: Currently compute_gfactors is only available when recon is
    %   is defined with recon = GenericCartesianGRAPPAReconstruction() and 
    %   not when recon = ImagesReconstructor() 
    
    % 3) set the reconstruction input to be the data we just preprocessed.
    recon.set_input(preprocessed_AcCont);
    
    % 4) Run the reconstruction using 'process' to call gadgetron.
    fprintf('---\n reconstructing...\n');
    recon.process();
    
    % Output
    
    % Reconstructed data sits in memory. We need to first get containers 
    % for the reconstructed images and g-factors, before extracting the 
    % data as MATLAB arrays. Containers in effect point to the data.
    
    % Get images and gfactors as containers with type mGadgetron.ImagesContainer
    images_imcont = recon.get_output('image');
    gfacts_imcont = recon.get_output('gfactor');
    
    % Return as MATLAB matrices the data pointed to by the containers.
    % Note the image data is complex.
    idata = images_imcont.as_array();
    gdata = gfacts_imcont.as_array();
    
    sl = 5 ; % Number of the slice to be displayed.
    if size(idata,3) < sl
        sl = 1 ;
    end
    
    % Display the modulus and phase for this reconstructed slice and gfactor.
    figure('Name',['Image and gfactor data, slice: ',num2str(sl)])
    subplot(2,2,1), imshow(abs(idata(:,:,sl)),[]), title('Image modulus')
    subplot(2,2,2), imshow(angle(idata(:,:,sl)),[-pi pi]), title('Image phase')
    
    subplot(2,2,3), imshow(abs(gdata(:,:,sl)),[]), title('Gfactor modulus')
    subplot(2,2,4), imshow(angle(gdata(:,:,sl)),[-pi pi]), title('Gfactor Phase')
    

%catch err
%    % display error information
%    fprintf('%s\n', err.message)
%    fprintf('error id is %s\n', err.identifier)
%end
