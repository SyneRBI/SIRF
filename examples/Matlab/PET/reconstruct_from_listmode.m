function reconstruct_from_listmode(engine)
% A demo showing reconstruction from raw data in listmode format.

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
import_str = set_up_PET(engine);
eval(import_str)

AcquisitionData.set_storage_scheme('memory');

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    MessageRedirector('info.txt', 'warn.txt');

    % create listmode-to-sinograms converter object
    lm2sino = ListmodeToSinograms();

    [filename, pathname] = uigetfile...
        ('*.l.hdr*', 'Select listmode data file', pet_data_path);
    list_file = fullfile(pathname, filename);
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', pet_data_path);
    tmpl_file = fullfile(pathname, filename);
    [filename, pathname] = uigetfile...
        ('*.n.hdr*', 'Select ECAT8 normalization file', pet_data_path);
    norm_file = fullfile(pathname, filename);
    
    % set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix('sinograms')
    lm2sino.set_template(tmpl_file)

    % set interval
    lm2sino.set_time_interval(0, 100)

    % set flags
    lm2sino.flag_on('store_prompts')
    lm2sino.flag_off('interactive')

    % set up the converter
    lm2sino.set_up()

    % convert
    lm2sino.process()

    % get access to the sinograms
    acq_data = lm2sino.get_output();
    % copy the acquisition data into a Python array
    acq_array = acq_data.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions: %d x %d x %d\n', acq_dim)
    z = uint16(acq_dim(3)/2);
    mUtilities.show_2D_array(acq_array(:,:,z), ...
        'acquisition data', 'tang. pos.', 'views');

    % create initial image estimate of dimensions and voxel sizes
    % compatible with the scanner geometry (included in the AcquisitionData
    % object acq_data) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0);
    image_array = image.as_array();
    fprintf('image dimensions: %d x %d x %d\n', size(image_array))

    % create acquisition sensitivity model from ECAT8 normalization data
    asm = AcquisitionSensitivityModel(norm_file);
    asm.set_up(acq_data);

    % select acquisition model that implements the geometric
    % forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix();
    acq_model.set_normalization(asm)
    acq_model.set_up(acq_data, image)

    % define objective function to be maximized as
    % Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_max_segment_num_to_process(1)

    % select Ordered Subsets Maximum A-Posteriori One Step Late as the
    % reconstruction algorithm (since we are not using a penalty, or prior, in
    % this example, we actually run OSEM);
    % this algorithm does not converge to the maximum of the objective function
    % but is used in practice to speed-up calculations
    num_subsets = 2;
    recon = OSMAPOSLReconstructor();
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_input(acq_data)

    % set up the reconstructor based on a sample image
    % (checks the validity of parameters, sets up objective function
    % and other objects involved in the reconstruction, which involves
    % computing/reading sensitivity image etc etc.)
    fprintf('setting up, please wait, may take a while...\n')
    recon.set_up(image)

    % set the initial image estimate
    recon.set_current_estimate(image)

    % in order to see the reconstructed image evolution
    % open up the user's access to the iterative process
    % rather than allow recon.reconstruct to do all job at once
    z = 20;
    num_subiterations = 1;
    for iter = 1 : num_subiterations
        fprintf('\n--------------------- Subiteration %d\n', iter)
        % perform an iteration
        recon.update_current_estimate()
        % display the current image
        image_array = recon.get_current_estimate().as_array();
        the_title = sprintf('iteration %d', iter);
        mUtilities.show_2D_array(image_array(:,:,z), the_title, 'x', 'y');
    end

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end


