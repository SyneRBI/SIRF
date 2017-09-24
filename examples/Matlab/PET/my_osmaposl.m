function image = my_osmaposl...
    (init_image, obj_fun, prior, filter, num_subsets, num_subiterations)
% Function implementing OSMAPOSL, an Ordered Subsets (OS) version of the One 
% Step Late algorithm (OSL) from Green et al for Maximum a Posteriori (MAP) 
% maximisation. 

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

if nargin < 1 % for run_all.m to ignore this function
    image = [];
    return
end

image = init_image;
for iter = 1 : num_subiterations
    fprintf('\n------------- Subiteration %d\n', iter) 

    % select subset
    subset = iter - 1;

    % get sensitivity as ImageData
    sens_image = obj_fun.get_subset_sensitivity(subset);

    % get backprojection of the ratio of measured to estimated
    % acquisition data
    grad_image = obj_fun.get_backprojection_of_acquisition_ratio...
                 (image, subset);

    % get gradient of prior as ImageData
    prior_grad_image = prior.get_gradient(image);

    % copy to Python arrays
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

    % fill current image with new values
    image.fill(image_array);

    % apply filter
    filter.apply(image);

end
end