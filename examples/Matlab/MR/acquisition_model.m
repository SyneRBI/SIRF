function acquisition_model(engine)
% USING_ACQUISITION_MODEL Demo illustrating creating simulated acquisitions 
% by the MR acquisition model's forward projection and their backprojecting
% to the image space
%
% In MATLAB, there are also ISMRMRD tools available for examining 
% data before processing.
%
% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
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

% Select and import SIRF MATLAB MR package so that SIRF MR objects can be 
% created in this function without using the prefix 'MR.'
if nargin < 1
    engine = [];
end
MR = set_up_MR(engine);
mr_data_path = sirf.Utilities.examples_data_path('MR');

AD = MR.AcquisitionData();

% acquisitions will be read from an HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = MR.AcquisitionData(fullfile(pathname, filename));

% pre-process acquisition data
fprintf('processing acquisitions...\n')
processed_data = MR.preprocess_acquisition_data(acq_data);
processed_data.sort()

% perform reconstruction to obtain a meaningful ImageData object
% (cannot be obtained in any other way at present):
% 1. Create a reconstruction object using 2D inverse Fourier transform and
if acq_data.is_undersampled()
%    CartesianGRAPPAReconstructor() sets up a default gadget chain.
    recon = MR.CartesianGRAPPAReconstructor();
    recon.compute_gfactors(false)
else
%    FullySampledCartesianReconstructor() sets up a default gadget chain.
    recon = MR.FullySampledCartesianReconstructor();
end
% 2. Provide pre-processed k-space data as input
recon.set_input(processed_data)

% 3. Run (i.e. 'process') the reconstruction.
fprintf('reconstructing...\n')
recon.process()

% retrieve reconstruction as ImageData object
image_data = recon.get_output();

par = {'flags', 'channels', 'slice', 'repetition'};
for i = 1 : image_data.number()
    fprintf('--- image %d\n', i)
    image = image_data.image(i);
    for p = 1 : length(par)
        name = char(par(p));
        fprintf('%s %d\n', name, image.info(name))
    end
    fprintf('matrix size %dx%dx%d\n', image.matrix_size())
end

% compute coil sensitivity maps
csms = MR.CoilSensitivityData();
csms.calculate(processed_data);

% create acquisition model based on the acquisition parameters
% stored in processed_data and image parameters stored in complex_images
acq_model = MR.AcquisitionModel(processed_data, image_data);
acq_model.set_coil_sensitivity_maps(csms);

% use the acquisition model (forward projection) to produce simulated
% acquisition data
simulated_acq_data = acq_model.forward(image_data);

% % test i/o
% simulated_acq_data.write('sim_acq_data.h5');
% simulated_acq_data = MR.AcquisitionData('sim_acq_data.h5');

% get simulated acquisition data as a Matlab double array
simulated_acq_array = simulated_acq_data.as_array();
% display simulated acquisition data
simulated_acq_array = permute(simulated_acq_array,[1,3,2]);
title = 'Simulated acquisition data (magnitude)';
sirf.Utilities.show_3D_array...
    (abs(simulated_acq_array).^0.2, title, 'samples', 'readouts', 'coil');
sirf.Utilities.set_window(0.1, 0.1, 0.8, 0.8)

% backproject simulated acquisition data
backprojected_data = acq_model.backward(simulated_acq_data);
% show backprojected data
backprojected_array = backprojected_data.as_array();
title = 'Backprojected data (magnitude)';
sirf.Utilities.show_3D_array...
    (abs(backprojected_array), title, 'samples', 'readouts', 'slice');

end
