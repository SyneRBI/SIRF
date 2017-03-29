function steepest_ascent(engine)
% Steepest ascent demo.
% Applies few steps of steepest ascent for the minimization of Poisson 
% logarithmic likelihood objective function using gradient for subset 0.

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
    % create uniform image
    image = ImageData();
    image_size = [111, 111, 31];
    voxel_size = [3, 3, 3.375];
    image.initialise(image_size, voxel_size)
    image.fill(1.0)

    % create filter that zeroes the image outside a cylinder of the same
    % diameter as the image xy-section size
    filter = TruncateToCylinderProcessor();
    % apply the filter to the image (zero the image outside the cylindric FOV)
    filter.apply(image)

    % define acquisition data
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file', pet_data_path);
    acq_data = AcquisitionData(fullfile(pathname, filename));

    % define the acquisition model
    acq_model = AcquisitionModelUsingRayTracingMatrix();

    % create objective function of Poisson logarithmic likelihood type
    % compatible with the acquisition data type
    obj_fun = make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_acquisition_data(acq_data)
    obj_fun.set_num_subsets(12)
    obj_fun.set_up(image)

    z = 20;
    
    % show the initial image
    image_array = image.as_array();
    figure(1)
    imshow(image_array(:,:,z));
    
    eps = 1.0e-6; % single precision round-off error level
    tau = 0.3; % steepest ascent step size
    
    v = obj_fun.value(image);
    fprintf('initial objective function value: %e\n', v);
    
    for iter = 1 : 3

        % obtain the gradient for subset 0
        grad = obj_fun.gradient(image, 0);
        filter.apply(grad)
        % zero the gradient outside the cylindric FOV
        grad_array = grad.as_array();

        max_image = max(max(max(image_array)));
        max_grad = max(max(max(abs(grad_array))));
        delta = max_grad*eps;

        % find maximal steepest descent step parameter t in image + t*grad 
        % such that the new image remains positive
        % since image is non-negative, the step size is limited by negative
        % gradients: it should not exceed -image/grad = abs(image/grad) at
        % points where grad is negative, thus, maximal t is min(abs(image/grad))
        % taken over such points

        % avoid division by zero at the next step
        grad_array(abs(grad_array) < delta) = delta;
        ratio = abs(image_array./grad_array);
        grad_array(image_array < 0) = 0;

        % select points inside cylindric FOV at which the gradient is negative
        select = image_array > 0 & grad_array < 0;
        if any(any(any(select)))
            % take the minimum of abs(image/grad) over selected points
            maxstep = min(min(min(ratio(select))));
        else
            % no such points - use a plausible value based on tau and
            % 'average' image-to-gradient ratio
            maxstep = tau*max_image/max_grad;
        end
        
        % at some voxels image values may be close to zero and the gradient may
        % also be close to zero there; hence, both may become negative because
        % of round-off errors;
        % find such voxels and freeze the image values at them
        exclude = image_array <= 0 & grad_array < 0;
        grad_array(exclude) = 0;
        
        t = min(maxstep, tau*max_image/max_grad);

        % perform steepest descent step
        fprintf('step %d, max change in image %e\n', iter, t*max_grad);
        image_array = image_array + t*grad_array;
        image.fill(image_array);
        filter.apply(image);

        % show the current image estimate
        image_array = image.as_array();
        figure(iter + 1)
        imshow(image_array(:,:,z));
        
        % quit if the image got substantially negative values
        min_image = min(min(min(image_array)));
        if min_image < -eps
            fprintf('image minimum is negative: %e\n', min_image)
            break
        end

    end

    v = obj_fun.value(image);
    fprintf('attained objective function value: %e\n', v);
    
catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end