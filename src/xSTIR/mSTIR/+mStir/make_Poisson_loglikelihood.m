function obj_fun = make_Poisson_loglikelihood(acq_data, features)
% Creates objective function.
% make_Poisson_loglikelihood(acq_data, features) returns the objective 
% function that works with acquisition data of the same kind as the first 
% argument and has features specified by the second argument;
% acq_data: AcquisitionData object
% features: Matlab string

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

    obj_fun = mStir.PoissonLogLikelihoodWithLinearModelForMeanAndProjData();
    obj_fun.set_acquisition_data(acq_data);
end