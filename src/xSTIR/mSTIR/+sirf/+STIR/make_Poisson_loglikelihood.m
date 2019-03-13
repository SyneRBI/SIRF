function obj_fun = make_Poisson_loglikelihood(acq_data)
% Creates an objective function for data with Poisson noise.
% Usage: 
%     obj_func = make_Poisson_loglikelihood(acq_data);
% acq_data: AcquisitionData object
% returns an objective function appropriate for the 
% specified acquisition data.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
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

    sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
    obj_fun = sirf.STIR.PoissonLogLikelihoodWithLinearModelForMeanAndProjData();
    obj_fun.set_acquisition_data(acq_data);
end