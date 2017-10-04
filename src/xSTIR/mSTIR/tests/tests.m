function tests(engine)
% GRAPPA reconstruction with the steepest descent step
% to illustrate the use of Acquisition Model projections.

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

% Select and import SIRF MATLAB MR package so that SIRF MR objects can be 
% created in this function without using the prefix 'MR.'
if nargin < 1
    engine = [];
end
import_str = set_up_PET(engine);
eval(import_str)

failed = 0;
eps = 1e-4;
ntest = 0;

% define raw data source
filename = 'Utahscat600k_ca_seg4.hs';
pathname = pet_data_path();
acq_data = AcquisitionData(fullfile(pathname, filename));
s = acq_data.norm();
v = variance(acq_data);
fprintf('---\n acquisition data norm: %e, variance: %e\n', s, v)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 3099.322, s, 0, eps);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 5.444323, v, 0, eps);

image = ImageData();
image_size = [111, 111, 31];
n = prod(image_size);
voxel_size = [3, 3, 3.375];
image.initialise(image_size, voxel_size)
image.fill(1.0)
filter = TruncateToCylinderProcessor();
filter.apply(image)
s = image.norm();
v = variance(image);
fprintf('---\n filtered image norm: %e, variance: %e\n', s, v)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 541.678, s, 0, eps);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0.178068, v, 0, eps);

prior = QuadraticPrior();
prior.set_penalisation_factor(0.5)

matrix = RayTracingMatrix();
matrix.set_num_tangential_LORs(2)

am = AcquisitionModelUsingMatrix();
am.set_matrix(matrix)

num_subsets = 12;

obj_fun = make_Poisson_loglikelihood(acq_data);
obj_fun.set_acquisition_model(am)
obj_fun.set_num_subsets(num_subsets)
obj_fun.set_up(image)

subset = 0;

% get sensitivity as ImageData
sens_image = obj_fun.get_subset_sensitivity(subset);

% get backprojection of the ratio of measured to estimated
% acquisition data
grad_image = obj_fun.get_backprojection_of_acquisition_ratio...
             (image, subset);

% get gradient of prior as ImageData
prior_grad_image = prior.get_gradient(image);

% copy to Matlab arrays
image_array = image.as_array();
sens_array = sens_image.as_array();
grad_array = grad_image.as_array();
prior_grad_array = prior_grad_image.as_array();

% update image data
denom = sens_array + prior_grad_array./num_subsets;
eps = 1e-6*max(abs(denom(:)));
denom(denom < eps) = eps; % avoid division by zero
update = grad_array./denom;
image_array = image_array.*update;

s = norm(update(:));
v = var(update(:));
delta = max(update(:))*eps;
fprintf('---\n update norm: %e, variance: %e\n', s, v)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 2377.229, s, 0, eps);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 14.7674, v, delta, eps);

s = norm(image_array(:));
v = var(image_array(:));
delta = max(image_array(:))*eps;
fprintf('---\n updated image norm: %e, variance: %e\n', s, v)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 7.610105, s, 0, eps);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0.000052, v, delta, eps);

if failed == 0
    fprintf('all tests passed\n')
else
    fprintf('%d tests failed\n', failed)
end
end

function v = variance(x)
a = double(x.as_array());
v = var(a(:));
end

function failed = test_failed(ntest, expected, actual, abstol, reltol)
failed = abs(expected - actual) > abstol + reltol*expected;
    if failed
        fprintf('+++ test %d failed\n', ntest)
    else
        fprintf('+++ test %d passed\n', ntest)
    end
end
