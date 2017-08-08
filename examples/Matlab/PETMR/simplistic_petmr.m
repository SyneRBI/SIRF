function simplistic_petmr(mr_engine, pet_engine)
% SIMPLISTIC_PETMR  Reconstruct fully sampled MR data then apply PET filter
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file in the ISMRMRD format.

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

% default engines to be used if none given
if nargin < 1
    mr_engine = [];
    pet_engine = [];
end
% emulate Python's import <module> as <alias>
import_mr_str = set_up_MR(mr_engine, 'MR');
eval(import_mr_str)
import_pet_str = set_up_PET(pet_engine, 'PET');
eval(import_pet_str)
% note that this will create copies of the engines module folders
% named +MR and +PET, so you must have sufficient privileges to write 
% to installation folder

% acquisitions will be read from an HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
input_data = MR.AcquisitionData(fullfile(pathname, filename));

% pre-process acquisition data
fprintf('processing acquisitions...\n')
processed_data = MR.preprocess_acquisition_data(input_data);

% perform MR reconstruction:
% 1. Create a reconstruction object using 2D inverse Fourier transform and
%    FullySampledCartesianReconstructor() sets up a default gadget chain.
recon = MR.FullySampledCartesianReconstructor();
% 2. Provide pre-processed k-space data as input
recon.set_input(processed_data)
% 3. Run (i.e. 'process') the reconstruction.
fprintf('reconstructing...\n')
recon.process()
% retrieve reconstruction as an ImageData object
mr_image = recon.get_output();

% display MR image
mUtilities.show_3D_array(abs(mr_image.as_array()), 'MR image data', ...
    'x (FE)', 'y (PE)', 'slice')

% convert MR image to PET image
image_array = abs(mr_image.as_array());
pet_image = PET.ImageData();
pet_image.initialise(size(image_array));
pet_image.fill(image_array);

% apply PET filter to PET image
filter = PET.TruncateToCylinderProcessor();
% standard data processor usage
filter.set_input(pet_image)
filter.process();
pet_image = filter.get_output();
% shortcuts for the above 3 lines
%pet_image = filter.process(pet_image);
%filter.apply(pet_image);

% display filtered PET image
mUtilities.show_3D_array(pet_image.as_array(), 'PET image data', 'x', 'y', 'slice')

end