function fully_sampled_tests(engine)
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
import_str = set_up_MR(engine);
eval(import_str)

failed = 0;
eps = 1e-4;
ntest = 0;

% define raw data source
filename = 'simulated_MR_2D_cartesian.h5';
pathname = mr_data_path();
acq_data = AcquisitionData(fullfile(pathname, filename));
acq_data_norm = acq_data.norm();
fprintf('---\n acquisition data norm: %e\n', acq_data_norm)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 221.201, acq_data_norm, 0, eps);

prep_gadgets = {'RemoveROOversamplingGadget'};
processed_data = acq_data.process(prep_gadgets);
processed_data_norm = processed_data.norm();
fprintf('---\n processed acquisition data norm: %e\n', processed_data_norm);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 209.021, processed_data_norm, 0, eps);

recon = FullySampledCartesianReconstructor();
recon.set_input(processed_data);
recon.process();
complex_images = recon.get_output();

complex_images_norm = complex_images.norm();
fprintf('---\n reconstructed images norm: %e\n', complex_images_norm);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 209.021, processed_data_norm, 0, eps);

csms = CoilSensitivityData();

fprintf('---\n sorting acquisitions...')
processed_data.sort()
fprintf('ok\n')
fprintf('---\n computing sensitivity maps...')
csms.calculate(processed_data)
fprintf('ok\n')

am = AcquisitionModel(processed_data, complex_images);

am.set_coil_sensitivity_maps(csms)

fwd_acqs = am.forward(complex_images);

fwd_acqs_norm = fwd_acqs.norm();
fprintf('---\n reconstructed images forward projection norm %e\n', ...
      fwd_acqs_norm)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 209.021, fwd_acqs_norm, 0, eps);

acqs_diff = fwd_acqs - processed_data;
rr = acqs_diff.norm()/fwd_acqs_norm;
fprintf('---\n reconstruction residual norm (rel): %e\n', rr)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0, rr, 1e-5, 0);

bwd_images = am.backward(processed_data);
imgs_diff = bwd_images - complex_images;
rd = imgs_diff.norm()/complex_images.norm();
fprintf(...
    '---\n difference between reconstructed and back-projected images: %e\n', ...
    rd)
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0, rd, 1e-5, 0);

xFy = processed_data * fwd_acqs;
Bxy = bwd_images * complex_images;
fprintf('---\n (x, F y) = (%e, %e)\n', real(xFy), imag(xFy))
fprintf('= (B x, y) = (%e, %e)\n', real(Bxy), imag(Bxy))
ntest = ntest + 1;
failed = failed + test_failed(ntest, real(xFy), real(Bxy), 0, 1e-5);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0, imag(xFy)/real(xFy), 1e-5, 0);
ntest = ntest + 1;
failed = failed + test_failed(ntest, 0, imag(Bxy)/real(Bxy), 1e-5, 0);

if failed == 0
    fprintf('all tests passed\n')
else
    fprintf('%d tests failed\n', failed)
end
end

function failed = test_failed(ntest, expected, actual, abstol, reltol)
failed = abs(expected - actual) > abstol + reltol*expected;
    if failed
        fprintf('+++ test %d failed\n', ntest)
    else
        fprintf('+++ test %d passed\n', ntest)
    end
end
