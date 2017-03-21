function fully_sampled_recon_three_chains
% FULLY_SAMPLED_RECON_THREE_CHAINS
% Runs 3 gadget chains of different type:
% - acquisition processing chain,
% - reconstruction chain,
% - image processing chain
% and how to visualise or modify data in between these chains.
%
% See also FULLY_SAMP_BASIC

% load mutilities and mgadgetron libraries
ccp_libload

% import mGadgetron MATLAB classes so that they can be called in this
% function without using the prefix 'mGadgetron.'
import mGadgetron.*

% acquisitions will be read from this HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
input_data = AcquisitionData(fullfile(pathname, filename));

% process data using Acquisitions processing chain
acq_proc = AcquisitionDataProcessor({'RemoveROOversamplingGadget'});
fprintf('processing acquisitions...\n')
preprocessed_data = acq_proc.process(input_data);

% As an example, here we access the preprocessed k-space and apply a
% Gaussian^4 weighting that will blur the image.
% Provides an example of interacting with the data between calls to Gadgetron.


data_array = preprocessed_data.as_array() ;
[nx, nc, ns] = size(data_array) ;
x = -(nx - 1)/2 : (nx - 1)/2;
sigma = 40;
w = exp(-x.*x/(2*sigma*sigma));
%w = window(@gausswin, nx) ;
w = reshape(w,[ nx 1 1]) ;
w = w.^4 ;
data_array = data_array .* repmat(w,[1 nc ns]) ;

% Re-fill the preprocessed_data object
preprocessed_data.fill(data_array)

% build reconstruction chain, here using a pre-set Set of gadgets. Can
% alternatively pass in list of gadgets as a cell array of gadget names.
recon = Reconstructor({'SimpleReconGadgetSet'});

% provide pre-processed k-space data to recon
recon.set_input(preprocessed_data)

% perform reconstruction
fprintf('reconstructing...\n')
recon.process()

% get reconstructed image object
complex_images = recon.get_output();

% extract real images using Images processing chain
% Note this still returns an mGadgetron.ImageData object that requires use
% of as_array() or show() to visulaise.
img_proc = ImageDataProcessor({'ExtractGadget'});
fprintf('processing images...\n')
images = img_proc.process(complex_images);

% plot obtained images
% See other demos for use of as_array() to extract a MATLAB array and then
% plot
images.show()

