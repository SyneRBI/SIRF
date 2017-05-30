function using_acquisition_model(engine)
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
import_str = setup_PET(engine);
eval(import_str)

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    MessageRedirector('info.txt');

    % create uniform image
    image = ImageData();
    image_size = [111, 111, 31];
    voxel_size = [3, 3, 3.375];
    image.initialise(image_size, voxel_size)

    % create a shape
    shape = EllipticCylinder();

    % add a shape to the image
    shape.set_length(400);
    shape.set_radii([40, 100]);
    shape.set_origin([60, 0, 10]);
    image.add_shape(shape, 1.0)

    % add another shape
    shape.set_radii([30, 30])
    shape.set_origin([-30, 60, 10])
    image.add_shape(shape, 1.5)

    % add another shape
    shape.set_origin([-30, -60, 10])
    image.add_shape(shape, 0.75)

    % z-coordinate of the xy-section to display
    z = int32(image_size(3)/2);

    % display the created phantom image
    image_array = image.as_array();
    mUtil.show_2D_array(image_array(:,:,z), 'phantom', 'x', 'y');

    % select the acquisition model that implements the geometric
    % forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix();
    % project the image to obtain simulated acquisition data;
    % raw data selected by the user is used as a template
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', pet_data_path);
    template = AcquisitionData(fullfile(pathname, filename));
    fprintf('setting up acquisition model...\n')
    acq_model.set_up(template, image)
    fprintf('projecting...\n')
    simulated_data = acq_model.forward(image); % simulated_data stored in memory
%     simulated_data = acq_model.forward(image, 'demo4data.hs'); % simulated_data stored in a file

    % display simulated data
    acq_array = simulated_data.as_array();
    acq_dim = size(acq_array);
    mUtil.show_2D_array(acq_array(:,:,acq_dim(3)/2), ...
        'simulated acquisition data', 'tang. pos.', 'views');

    % backproject the simulated data
    fprintf('backprojecting...\n')
    backprojected_image = acq_model.backward(simulated_data);
    % display backprojected data
    image_array = backprojected_image.as_array();
    mUtil.show_2D_array(image_array(:,:,z), ...
        'backprojection of simulated data', 'x', 'y');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end