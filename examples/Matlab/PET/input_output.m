function input_output(engine)
% Acquisition model demo: creates an image, forward projects it to simulate
% acquisition data and then backprojects the result.

% Uncomment the next line to get more information on PET acquisition model
% doc AcquisitionModel

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
% import_str = set_up_PET(engine);
% eval(import_str)
PET = set_up_PET(engine);

try
    % create acquisition data from scanner parameters
    fprintf('creating acquisition data...\n')
    acq_data = PET.AcquisitionData('Siemens_mMR');
    % set all values to 1.0
    acq_data.fill(1.0);

    % copy the acquisition data into a Python array and display
    acq_array = acq_data.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions: %d x %d x %d\n', acq_dim)
    z = uint16(acq_dim(3)/2);
    sirf.Utilities.show_2D_array(acq_array(:,:,z), 'Acquisition data',...
        'tang. pos.', 'views');

    % create image of dimensions and voxel sizes compatible with the scanner
    % geometry (stored in the AcquisitionData object ad)
    % and initialize each voxel to 2.0
    image = acq_data.create_uniform_image(2.0);
    % show the image
    image_array = image.as_array();
    image_dim = size(image_array);
    fprintf('image dimensions: %d x %d x %d\n', image_dim)
    z = uint16(image_dim(3)/2);
    sirf.Utilities.show_2D_array(image_array(:,:,z), 'Image', 'x', 'y');

    % write acquisition data and image to files
    fprintf('writing acquisition data...\n')
    acq_data.write('ones');
    fprintf('writing image...\n')
    image.write('twos');

    % read acquisition data and image from files and display
    acq = PET.AcquisitionData('ones.hs');
    acq_array = acq.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions: %d x %d x %d\n', acq_dim)
    z = uint16(acq_dim(3)/2);
    sirf.Utilities.show_2D_array(acq_array(:,:,z), 'Acquisition data',...
        'tang. pos.', 'views');
    img = PET.ImageData();
    img.read_from_file('twos.hv');
    image_array = img.as_array();
    image_dim = size(image_array);
    fprintf('image dimensions: %d x %d x %d\n', image_dim)
    z = uint16(image_dim(3)/2);
    sirf.Utilities.show_2D_array(image_array(:,:,z), 'Image', 'x', 'y');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end