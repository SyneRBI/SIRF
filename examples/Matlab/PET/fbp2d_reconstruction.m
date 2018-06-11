function fbp2d_reconstruction(engine)
% FBP2D reconstruction demo. 

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

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    MessageRedirector('info.txt');

    % PET acquisition data to be read from this file
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file', pet_data_path);
    acq_data = AcquisitionData(fullfile(pathname, filename));
    
    % create reconstructor object
    recon = FBP2DReconstructor();
    % specify the acquisition data

    % reconstruct with default settings
    recon.set_input(acq_data);
    recon.reconstruct();
    image = recon.get_output();
    dim = image.size();
    fprintf('image size: %d x %d x %d\n', dim(1), dim(2), dim(3));
    z = idivide(2*dim(3), 3);
    % display the reconstructed image
    image.show(z);
    
    % change image size
    recon.set_output_image_size_xy(2*dim(1));
    recon.reconstruct();
    image = recon.get_output();
    image.show(z);

    % zoom in
    zoom = 2.5;
    recon.set_zoom(zoom);
    recon.reconstruct();
    image = recon.get_output();
    image.show(z);

    % use a Hann filter
    alpha = 0.5;
    recon.set_alpha_cosine_window(alpha);
    recon.reconstruct();
    image = recon.get_output();
    image.show(z);

    % a Hann filter with lower cut-off (0.5 is no cut-off)
    fc = 0.2;
    recon.set_frequency_cut_off(fc);
    recon.reconstruct();
    image = recon.get_output();
    image.show(z);

    % alternative way to set the output image parameters (via image template)
    image_tmpl = acq_data.create_uniform_image(); % image template
    recon.set_up(image_tmpl); % use image template to create the output image
    recon.reconstruct();
    image = recon.get_output();
    image.show(z);

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end