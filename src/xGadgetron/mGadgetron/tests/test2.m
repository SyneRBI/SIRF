function [failed, ntests] = test2(record, engine)
% MR test set 2.

% GRAPPA reconstruction and the illustration of the use of 
% Acquisition Model projections.

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
if nargin < 2
    engine = [];
end
if nargin < 1
    record = false;
end
% import_str = set_up_MR(engine);
% eval(import_str)
MR = set_up_MR(engine);

test = mUtilities.mTest('test2.txt', record);

filename = 'simulated_MR_2D_cartesian_Grappa2.h5';
pathname = mUtilities.examples_data_path('MR');
acq_data = MR.AcquisitionData(fullfile(pathname, filename));
test.check(acq_data.norm())

prep_gadgets = {'RemoveROOversamplingGadget'};
processed_data = acq_data.process(prep_gadgets);
test.check(processed_data.norm())

recon = MR.CartesianGRAPPAReconstructor();
recon.compute_gfactors(false);
recon.set_input(processed_data);
recon.process();
complex_images = recon.get_output();
test.check(complex_images.norm())

processed_data.sort()
csms = MR.CoilSensitivityData();
csms.calculate(processed_data)

am = MR.AcquisitionModel(processed_data, complex_images);
am.set_coil_sensitivity_maps(csms)
fwd_acqs = am.forward(complex_images);
fwd_acqs_norm = fwd_acqs.norm();
test.check(fwd_acqs_norm)

acqs_diff = fwd_acqs - processed_data;
rr = acqs_diff.norm()/fwd_acqs_norm;
test.check(rr, 1e-4)

bwd_images = am.backward(processed_data);
imgs_diff = bwd_images - complex_images;
rd = imgs_diff.norm()/complex_images.norm();
test.check(rd, 1e-4)

xFy = processed_data * fwd_acqs;
Bxy = bwd_images * complex_images;
test.check(abs(real(xFy)/real(Bxy) - 1.0), 1e-4)

failed = test.failed;
ntests = test.ntest;

if record
    fprintf('%d measurements recorded\n', ntests)
elseif failed == 0
    fprintf('all %d tests passed\n', ntests);
else
    fprintf('%d out of %d tests failed\n', failed, ntests);
end
end