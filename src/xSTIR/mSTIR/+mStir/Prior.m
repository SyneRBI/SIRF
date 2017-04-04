classdef Prior < handle
% Class for objects handling the prior: a penalty term to be added to the
% objective function maximized by iterative reconstruction algorithms.

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

    properties
        name
        handle
    end
    methods
        function self = Prior()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_penalisation_factor(self, value)
%***SIRF*** Sets the factor by which the penalty term (prior) is to be multiplied
%         before adding to the objective function.
            mStir.setParameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', value, 'f')
        end
        function value = get_penalisation_factor(self)
%***SIRF*** Returns the penalty factor in front of the prior.
            value = mStir.parameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', 'f');
        end
        function grad = get_gradient(self, image)
%***SIRF*** Returns the value of the gradient of the prior for the specified 
%         image.
            grad = mStir.ImageData();
            grad.handle = calllib('mstir', 'mSTIR_priorGradient', ...
                self.handle, image.handle);
            mUtil.checkExecutionStatus('Prior:get_gradient', grad.handle)
        end
    end
end