function fully_sampled_recon(engine)
% FULLY_SAMPLED_RECON  Recon of fully sampled data
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file in the ISMRMRD format.
%    Example datasets:
%    a) 'meas_MID00103_FID57244_test.dat' is 
%       available from https://www.ccpsynerbi.ac.uk/downloads
%       This is in the manufacturer's raw data format and needs to be
%       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
%       This executable is installed on the Virtual Machine.
%       This is a large dataset from a fully sampled, 3D acquisition.
% 
%    b) 'meas_MID00107_FID57248_test_2D.dat' similar to above but is
%    multi-slce 2D (currently not avaialable for download). Slices not
%    acquired in sequential order.
%
%    c) An ISMRMRD h5 file can be simulated using test_create_dataset.m
%    available in the ISMRMRD/examples/matlab folder at:
%    https://github.com/ismrmrd/ismrmrd
%    This gives 5 repetitions of a single slice.
%
% See also FULLY_SAMPLED_RECON_SINGLE_CHAIN
% FULLY_SAMPLED_RECON_THREE_CHAINS

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2019 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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

% MR raw data formats from different vendors can be transformed to 
% HDF file format using siemens_to_ismrmrd, philips_to_ismrmrd or
% bruker_to_ismrmrd on https://github.com/ismrmrd/.
% acquisitions will be read from this HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
input_data = MR.AcquisitionData(fullfile(pathname, filename));

% pre-process acquisition data
% Prior to image reconstruction several pre-processing steps such as 
% asymetric echo compensation, noise decorrelation for multi-coil data or 
% removal of oversampling along frequency encoding (i.e. readout or kx)
% direction. So far only the removal of readout oversampling and noise and
% asymmetric echo adjusting is implemented
fprintf('processing acquisitions...\n')
processed_data = MR.preprocess_acquisition_data(input_data);

% perform reconstruction:
% 1. Create a reconstruction object using 2D inverse Fourier transform and
%    FullySampledCartesianReconstructor() sets up a default gadget chain.
recon = MR.FullySampledCartesianReconstructor();
% 2. Provide pre-processed k-space data as input
recon.set_input(processed_data)

% 3. Run (i.e. 'process') the reconstruction.
fprintf('reconstructing...\n')
recon.process()

% retrieve reconstruction as ImageDataobject
image_data = recon.get_output();

% show reconstructed image data
if exist('montage','file') && exist('mat2gray','file')
    % use as_array method to copy image_data to a Matlab array
    idisp = mat2gray(abs(image_data.as_array()));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
else
    sirf.Utilities.show_3D_array...
        (abs(image_data.as_array()), 'Reconstructed image data (magnitude)', ...
        'samples', 'readouts', 'slice');
end

% filter image
image_array = image_data.as_array();
image_array(abs(image_array) < 0.2*max(image_array(:))) = 0;
image_data.fill(image_array);
sirf.Utilities.show_3D_array...
    (abs(image_data.as_array()), 'Filtered image data (magnitude)', ...
    'samples', 'readouts', 'slice');
