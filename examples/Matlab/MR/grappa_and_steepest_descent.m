function grappa_and_steepest_descent(engine)
% GRAPPA reconstruction with the steepest descent step
% to illustrate the use of Acquisition Model projections.

% Select and import SIRF MATLAB MR package so that SIRF MR objects can be 
% created in this function without using the prefix 'MR.'
% set_up_mr
% import MR.*
if nargin < 1
    engine = [];
end
eval(setup_MR(engine))

% define raw data source
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
input_data = AcquisitionData(fullfile(pathname, filename));

% pre-process acquisitions
fprintf('---\n preprocessing...\n');
preprocessed_data = preprocess_acquisitions(input_data);
pp_norm = preprocessed_data.norm();

% perform reconstruction
recon = GenericCartesianGRAPPAReconstruction();
recon.compute_gfactors(false);
recon.set_input(preprocessed_data);
fprintf('---\n reconstructing...\n');
recon.process();
complex_images = recon.get_output();

% compute coil sensitivity maps
csms = CoilSensitivityMaps();
fprintf('---\n sorting acquisitions...\n')
preprocessed_data.sort()
fprintf('---\n calculating sensitivity maps...\n')
csms.calculate(preprocessed_data)

% create acquisition model based on the acquisition parameters
% stored in preprocessed_data and image parameters stored in complex_images
am = AcquisitionModel(preprocessed_data, complex_images);
am.set_coil_sensitivity_maps(csms)

% use the acquisition model (forward projection) to simulate acquisitions
fwd_data = am.forward(complex_images);
fwd_norm = fwd_data.norm();
% compute the difference between real and simulated acquisitions
diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm);
rr = diff.norm()/fwd_norm;
fprintf('---\n reconstruction residual norm (rel): %e\n', rr)

% try to improve the reconstruction by the steepest descent step
g = am.backward(diff);
w = am.forward(g);
alpha = (g*g)/(w*w);
r_complex_imgs = complex_images - g*alpha;

idata = abs(complex_images.as_array());
rdata = abs(r_complex_imgs.as_array());
idata = idata/max(max(max(idata)));
rdata = rdata/max(max(max(rdata)));
n = complex_images.number();

fprintf('---\nEnter slice number to view it\n')
fprintf('(a value outside the range [1 : %d] will stop this loop)\n', n)
while (true)
    i = input('slice: ');
    if i < 1 || i > n
        break
    end
    figure(i)
    imshow(idata(:,:,i));
    title('GRAPPA image')
    figure(i + n)
    imshow(rdata(:,:,i));
    title('refined image')
end


