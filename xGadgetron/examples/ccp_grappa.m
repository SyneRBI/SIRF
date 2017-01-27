function ccp_grappa
% CCP_GRAPPA Annotated demo for GRAPPA reconstruction
% 
% Demonstrates use of CCP SIRF code to perform a GRAPPA reconstruction.
%
% This version assumes the more modern gadget names available  
%  e.g. GenericReconFieldOfViewAdjustmentGadget 
%
% Usage:
%  Convert scanner raw data to ISMRMRD format e.g. by using the executable
%  siemens_to_ismrmrd. Example raw data is available from 
%  https://www.ccppetmr.ac.uk/downloads
%  Choose the GRAPPA dataset for this demo.
%
%  Ensure there is a listening Gadgetron - typically started in a terminal.
%
%  ccp_grappa
%
%
% Adapted by David Atkinson from original code by Evgueni Ovtchinnikov
% See also CCP_LIBLOAD


% load mutilities and mgadgetron libraries if not already loaded
ccp_libload

% Here choose Gadgetron as MR recon engine
import mGadgetron.*

try 
    % Get the filename of the input ISMRMRD h5 file
    filein = pref_uigetfile('ccp','filename');
    
    % Create mGadgetron.AcquisitionData
    input_data = AcquisitionData(filein);
    
    % Pre-process acquisition data using Gadgetron
    
    % List preparation gadgets
    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoAdjustROGadget'} ...
         {'RemoveROOversamplingGadget'}];
    
    % Call Gadgetron by using the 'process' method. This runs the gadgets 
    % specified in prep_gadgets, returning an instance 
    % of an mGadgetron.AcquisitionsContainer
    preprocessed_AcCont = input_data.process(prep_gadgets);
    
    
    % Perform reconstruction
    
    % Create a recon object for the desired reconstruction, here 
    % with type mGadgetron.GenericCartesianGRAPPAReconstruction 
    %  (To find what this does behind the scenes:
    %  Type edit mGadgetron.GenericCartesianGRAPPAReconstruction
    %  and look for the name assigned (here SimpleGRAPPAReconstructionProcessor)
    %  Then find the gadget chain defined by the class with the same 
    %  name in the file xGadgetron/cGadgetron/chain_lib.h )
    recon = GenericCartesianGRAPPAReconstruction();
    
    % For undersampled acquisition data, the GRAPPA gadget can 
    % compute G-factors in addition to reconstructed images.
    recon.compute_gfactors(true);
    
    % Set the input acquisition for the recon
    recon.set_input(preprocessed_AcCont);
    
    fprintf('---\n reconstructing...\n');
    
    % Call the recon (here calls gadgetron again)
    recon.process();
    
    % Get images and gfactors as containers with type mGadgetron.ImagesContainer
    images_imcont = recon.get_output('image');
    gfacts_imcont = recon.get_output('gfactor');
    
    % Return as MATLAB matrices the data pointed to by the containers.
    % Note the image data is complex.
    idata = images_imcont.as_array();
    gdata = gfacts_imcont.as_array();
    
    % Display using eshow if available
    if exist('eshow','file')
        eshow(idata, 'Name','Complex recon images')
        eshow(gdata,' Name', 'G-factors')
    else
        % use show method on images container
        images_imcont.show ;
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
