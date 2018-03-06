function reconstruct_from_listmode(engine)
% A demo showing reconstruction from raw data in listmode format.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 Rutherford Appleton Laboratory STFC.
% Copyright 2018 University College London.
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
    % direct all information printing to info.txt; warnings to warn.txt
    % error messages to go to Matlab Command Window
    MessageRedirector('info.txt', 'warn.txt');

    % First step is to create AcquisitionData ("sinograms") from the
    % listmode file.
    % See the listmode_to_sinograms demo for some more information on this
    % step.

    % create listmode-to-sinograms converter object
    lm2sino = ListmodeToSinograms();

    [filename, pathname] = uigetfile...
        ('*.l.hdr', 'Select listmode data file', pet_data_path);
    list_file = fullfile(pathname, filename);
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', pet_data_path);
    tmpl_file = fullfile(pathname, filename);
    [filename, pathname] = uigetfile...
        ('*.n.hdr', 'Select ECAT8 normalization file', pet_data_path);
    norm_file = fullfile(pathname, filename);
    [filename, pathname] = uigetfile...
        ('*.*hv', 'Select attenuation file', pet_data_path);
    attn_file = fullfile(pathname, filename);
    
    % set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix('sinograms')
    lm2sino.set_template(tmpl_file)

    % set interval
    lm2sino.set_time_interval(0, 50)

    % set up the converter
    lm2sino.set_up()

    % convert
    fprintf('converting raw data to sinograms...\n')
    lm2sino.process()

    % estimate randoms
    fprintf('estimating randoms...\n')
    randoms = lm2sino.estimate_randoms();

    % get access to the sinograms
    acq_data = lm2sino.get_output();
    % copy the acquisition data into a Matlab array
    acq_array = acq_data.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions: %d x %d x %d\n', acq_dim)
    z = round(acq_dim(3)/2);
    mUtilities.show_2D_array(acq_array(:,:,z), ...
        'acquisition data', 'tang. pos.', 'views');

    % read attenuation image
    attn_image = ImageData(attn_file);
    attn_image_as_array = attn_image.as_array();
    % select a slice appropriate for the NEMA acquistion data
    z = 72;
    % z = round(size(attn_image_as_array, 3)/2);
    mUtilities.show_2D_array(attn_image_as_array(:,:,z), ...
        'attenuation image', 'tang. pos.', 'views');

    % create initial image estimate of dimensions and voxel sizes
    % compatible with the scanner geometry (included in the AcquisitionData
    % object acq_data) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0, 127, 127);
    image_array = image.as_array();
    fprintf('image dimensions: %d x %d x %d\n', size(image_array))

    % select acquisition model that implements the geometric
    % forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix();
    acq_model.set_num_tangential_LORs(10)

    % create acquisition sensitivity model from ECAT8 normalization data
    asm_norm = AcquisitionSensitivityModel(norm_file);
    % create acquisition sensitivity model for attenuation
    asm_attn = AcquisitionSensitivityModel(attn_image, acq_model);
    asm_attn.set_up(acq_data);
    % compute attenuation factors
    bin_eff = AcquisitionData(acq_data);
    bin_eff.fill(1.0);
    fprintf('applying attenuation (please wait, may take a while)...\n')
    asm_attn.unnormalise(bin_eff);
    %store these in a new acquisition sensitivity model
    asm_beff = AcquisitionSensitivityModel(bin_eff);

    % chain attenuation and ECAT8 normalisation
    asm = AcquisitionSensitivityModel(asm_norm, asm_beff);

    acq_model.set_acquisition_sensitivity(asm);
    acq_model.set_background_term(randoms);

    % define objective function to be maximized as a
    % log of the Poisson likelihood
    obj_fun = make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    % reduce number of segments that we will handle to save some time in
    % the demo
    obj_fun.set_max_segment_num_to_process(1)

    % select Ordered Subsets Maximum A-Posteriori One Step Late as the
    % reconstruction algorithm (since we are not using a penalty, or prior, in
    % this example, we actually run OSEM);
    % this algorithm does not converge to the maximum of the objective function
    % but is used in practice to speed-up calculations
    % See the reconstruction demos for more complicated examples     
    num_subsets = 7;
    num_subiterations = 2;
    recon = OSMAPOSLReconstructor();
    recon.set_objective_function(obj_fun);
    recon.set_num_subsets(num_subsets);
    recon.set_num_subiterations(num_subiterations);

    % set up the reconstructor based on a sample image
    % (checks the validity of parameters, sets up objective function
    % and other objects involved in the reconstruction, which involves
    % computing/reading sensitivity image etc etc.)
    fprintf('setting up reconstructor, please wait, may take a while...\n');
    recon.set_up(image);

    % set the initial image estimate
    recon.set_current_estimate(image);

    % reconstruct
    fprintf('reconstructing, please wait...\n')
    recon.process();

    % display the reconstructed image
    image_array = recon.get_current_estimate().as_array();
    the_title = sprintf('Reconstructed image');
    mUtilities.show_2D_array(image_array(:,:,z), the_title, 'x', 'y');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end


