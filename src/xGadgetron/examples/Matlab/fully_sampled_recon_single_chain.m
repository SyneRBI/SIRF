function fully_sampled_recon_single_chain
% FULLY_SAMPLED_RECON_SINGLE_CHAIN Complete recon in one Gadgetron chain
%
% See fully_samp_basic for example data files.
%
% See also FULLY_SAMP_BASIC  FULLY_SAMPLED_RECON


% load the SIRF mutilities library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

% This demo uses gadgetron as the MR reconstruction engine. 
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

import mGadgetron.*


% In this demo, rather than using a predefined image reconstruction 
% object, here an image reconstruction object is created by concatinating 
% multiple gadgets - for more information on Gadgetron and its gadgets,
% please see: 
%  https://github.com/gadgetron/.

% Parameters for individual gadgets can be defined either during the 
% creation of the reconstruction object:
% e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
% or by giving a gadget a label (cf. label ex: for the last gadget)
% and using set_gadget_property(label, propery, value).
% The gadgets will be concatenated and will be executed as soon as 
% process() is called.
gadgets = [...
    {'RemoveROOversamplingGadget'}, ...
    {'AcquisitionAccumulateTriggerGadget'}, ...
    {'BucketToBufferGadget'}, ...
    {'SimpleReconGadget'}, ...
    {'ImageArraySplitGadget'}, ...
    {'ex:ExtractGadget'} ...
    ];

% create reconstructor
recon = ImagesReconstructor(gadgets);

% change a property of the gadget labelled 'ex'
% ExtractGadget defines which type of image should be returned:
% none      0
% magnitude 1
% real      2
% imag      4
% phase     8
% max       16  
% in this example '5' returns both magnitude and imag 
recon.set_gadget_property('ex', 'extract_mask', 5);

% define raw data source
[filename, pathname] = uigetfile('*.h5', 'Select raw data file');
input_data = AcquisitionData(fullfile(pathname, filename));
recon.set_input(input_data)

% perform reconstruction
recon.process()

% get reconstructed image object
images = recon.get_output();

% plot obtained images (note imag are zero)
data = images.as_array();
figure('Name','images')
idisp = mat2gray(abs(data));
montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;


