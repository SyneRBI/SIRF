function fully_sampled_recon_three_chains(engine)
% FULLY_SAMPLED_RECON_THREE_CHAINS
% Runs 3 gadget chains of different type:
% - acquisition processing chain,
% - reconstruction chain,
% - image processing chain
% and how to visualise or modify data in between these chains.
%
% See also FULLY_SAMPLED_RECON

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% default engine to be used if none given
if nargin < 1
    engine = [];
end
% import_str = set_up_MR(engine);
% eval(import_str)
MR = set_up_MR(engine);
mr_data_path = sirf.Utilities.examples_data_path('MR');

% acquisitions will be read from this HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = MR.AcquisitionData(fullfile(pathname, filename));

% process data using Acquisitions processing chain
acq_proc = MR.AcquisitionDataProcessor({'RemoveROOversamplingGadget'});
fprintf('processing acquisitions...\n')
acq_proc.set_input(acq_data)
acq_proc.process();
preprocessed_data = acq_proc.get_output();
% a shortcut for the above 3 lines
% preprocessed_data = acq_proc.process(acq_data);

% As an example, here we access the preprocessed k-space and apply a
% Gaussian^4 weighting that will blur the image.
% Provides an example of interacting with the data between calls to Gadgetron.


preprocessed_data_array = preprocessed_data.as_array() ;
[nx, nc, ns] = size(preprocessed_data_array) ;
x = -(nx - 1)/2 : (nx - 1)/2;
sigma = 40;
w = exp(-x.*x/(2*sigma*sigma));
%w = window(@gausswin, nx) ;
w = reshape(w,[ nx 1 1]) ;
w = w.^4 ;
preprocessed_data_array = preprocessed_data_array .* repmat(w,[1 nc ns]) ;

% Re-fill the preprocessed_data object
preprocessed_data.fill(preprocessed_data_array)

% build reconstruction chain, here using a pre-set Set of gadgets. Can
% alternatively pass in list of gadgets as a cell array of gadget names.
recon = MR.Reconstructor({'SimpleReconGadgetSet'});

% provide pre-processed k-space data to recon
recon.set_input(preprocessed_data)

% perform reconstruction
fprintf('reconstructing...\n')
recon.process()

% get reconstructed image object
complex_image_data = recon.get_output();

% extract real images using Images processing chain
% Note this still returns an sirf.Gadgetron.ImageData object that requires use
% of as_array() or show() to visulaise.
img_proc = MR.ImageDataProcessor({'ExtractGadget'});
fprintf('processing images...\n')
img_proc.set_input(complex_image_data)
img_proc.process();
real_image_data = img_proc.get_output();
% a shortcut for the above 3 lines
% real_image_data = img_proc.process(complex_image_data);

% show obtained images
% See other demos for use of as_array() to extract a MATLAB array and then
% plot
title = 'Reconstructed image data (magnitude)';
sirf.Utilities.show_3D_array(abs(real_image_data.as_array()), title, ...
    'samples', 'readouts', 'slice');

