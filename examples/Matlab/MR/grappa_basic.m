function grappa_basic(engine)
% GRAPPA_BASIC Demo for reconstruction of GRAPPA acquired data.
% 
% Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
% See function grappa_detail.m for an example showing more of the 
% workings and functionality of the SIRF code with the explicit use of the
% Gadgetron reconstruction engine.
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
% The argument:
% engine: Matlab string with the name of the reconstruction package to be 
%         used, defaults to Gadgetron
%
% Usage:
%  grappa_basic
%  grappa_basic(engine)
%
%
% See also GRAPPA_DETAIL GEN_US_DATA

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

if nargin < 1
    engine = [];
end
% import_str = set_up_MR(engine);
% eval(import_str)
MR = set_up_MR(engine);
mr_data_path = sirf.Utilities.examples_data_path('MR');

% Get the filename of the input ISMRMRD h5 file
[fn,pn] = uigetfile('*.h5','Select ISMRMRD H5 file', mr_data_path) ;
filein = fullfile(pn,fn) ;

% Load this ISMRMRD h5 file, creating an input Container
acq_data = MR.AcquisitionData(filein);

% Pre-process this input data. (Currently this is a MATLAB script that just
% sets up a 3 chain gadget. In the future it will be independent of the MR
% recon engine.)
preprocessed_data = MR.preprocess_acquisition_data(acq_data);

% Perform reconstruction of the preprocessed data.
% 1. set the reconstruction to be for Cartesian GRAPPA data.
recon = MR.CartesianGRAPPAReconstructor();

% 2. set the reconstruction input to be the data we just preprocessed.
recon.set_input(preprocessed_data);

% 3. run (i.e. 'process') the reconstruction.
fprintf('---\n reconstructing...\n');
recon.process();

% Extract an image Container from the reconstruction and convert this
% to a MATLAB array.
image_data = recon.get_output('image');
image_array = image_data.as_array();  % returns a complex array

sl = 5 ; % Number of the slice to be displayed.
if size(image_array,3) < sl
    sl = 1 ;
end

% Display the modulus and phase for this reconstructed slice.
figure('Name',['idata, slice: ',num2str(sl)])
subplot(1,2,1), imshow(abs(image_array(:,:,sl)),[]), title('Abs')
subplot(1,2,2), imshow(angle(image_array(:,:,sl)),[-pi pi]), title('Phase')