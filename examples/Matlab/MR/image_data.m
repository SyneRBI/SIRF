function image_data(engine)
% USING_ACQUISITION_DATA Demo illustrating acquisitions pre-processing 
% and displaying.
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

% read image from an HDF file
[filename, pathname] = uigetfile('*.h5', 'Select image data file', mr_data_path);
img_data = MR.ImageData(fullfile(pathname, filename));
[nx, ny, nz] = img_data.dimensions();
fprintf('image data dimensions: %dx%dx%d\n', nx, ny, nz)
img_data.show()
img_copy = img_data.clone();
diff = img_data - img_copy;
fprintf('|image - image_copy| = %f\n', diff.norm())

% create image from acquisitions
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = MR.AcquisitionData(fullfile(pathname, filename));
[ns, nc, na] = acq_data.dimensions();
fprintf('raw data dimensions: %d samples, %d coils, %d acquisitions\n', ...
    ns, nc, na)
img_data = MR.ImageData(acq_data);
[nx, ny, nz] = img_data.dimensions();
fprintf('image data dimensions: %dx%dx%d\n', nx, ny, nz)
