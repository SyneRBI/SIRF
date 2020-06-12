function [failed, ntests] = test1(record, engine)
% PET test set 1.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2020 University College London.
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
if nargin < 2
    engine = [];
end
if nargin < 1
    record = false;
end
% import_str = set_up_PET(engine);
% eval(import_str)
PET = set_up_PET(engine);

test = sirf.Utilities.mTest('test1.txt', record);

% define raw data source
filename = 'Utahscat600k_ca_seg4.hs';
pathname = sirf.Utilities.examples_data_path('PET');
acq_data = PET.AcquisitionData(fullfile(pathname, filename));
s = acq_data.norm();
v = variance(acq_data);
test.check(s)
test.check(v)
disp("Printing AcqData info")
disp(acq_data.get_info())

image = PET.ImageData();
image_size = [111, 111, 31];
voxel_size = [3, 3, 3.375];
image.initialise(image_size, voxel_size)
image.fill(1.0)
filter = PET.TruncateToCylinderProcessor();
filter.apply(image)
s = image.norm();
v = variance(image);
test.check(s)
test.check(v)

prior = PET.QuadraticPrior();
prior.set_penalisation_factor(0.5)
prior.set_up(image);

matrix = PET.RayTracingMatrix();
matrix.set_num_tangential_LORs(2)

am = PET.AcquisitionModelUsingMatrix();
am.set_matrix(matrix)

num_subsets = 12;

obj_fun = PET.make_Poisson_loglikelihood(acq_data);
obj_fun.set_acquisition_model(am)
obj_fun.set_num_subsets(num_subsets)
obj_fun.set_up(image)

subset = 0;

sens_image = obj_fun.get_subset_sensitivity(subset);

grad_image = obj_fun.get_backprojection_of_acquisition_ratio...
             (image, subset);

prior_grad_image = prior.get_gradient(image);

image_array = image.as_array();
sens_array = sens_image.as_array();
grad_array = grad_image.as_array();
prior_grad_array = prior_grad_image.as_array();

denom = sens_array + prior_grad_array./num_subsets;
delta = 1e-6*max(abs(denom(:)));
denom(denom < delta) = delta; % avoid division by zero
update = grad_array./denom;
image_array = image_array.*update;

eps = 1e-4;
s = norm(update(:));
v = var(update(:));
delta = max(update(:))*eps;
test.check(s)
test.check(v, delta)

s = norm(image_array(:));
v = var(image_array(:));
delta = max(image_array(:))*eps;
test.check(s)
test.check(v, delta)

% Check that the storage scheme for a pre-initialised object
% doesn't change just because the default has changed.
temp = PET.AcquisitionData();
temp.set_storage_scheme('memory');
ad = PET.AcquisitionData(fullfile(pathname, filename));
temp.set_storage_scheme('file');
test.check_if_equal('memory', ad.get_storage_scheme());

failed = test.failed;
ntests = test.ntest;

if record
    fprintf('%d measurements recorded\n', ntests)
elseif failed == 0
    fprintf('all %d tests passed\n', ntests)
else
    fprintf('%d out of %d tests failed\n', failed, ntests)
end
end

function v = variance(x)
a = double(x.as_array());
v = var(a(:));
end
