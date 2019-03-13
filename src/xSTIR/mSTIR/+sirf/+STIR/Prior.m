classdef Prior < handle
% Class for objects handling the prior.
% The prior is a penalty term added to the objective function maximized 
% by iterative reconstruction algorithms.

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

    properties
        name
        handle_
    end
    methods (Static)
        function name = class_name()
            name = 'Prior';
        end
    end
    methods
        function self = Prior()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
            end
        end
        function set_penalisation_factor(self, value)
%***SIRF*** Sets the factor by which the penalty term (prior) is to be 
%         multiplied before adding to the objective function.
            sirf.STIR.setParameter...
                (self.handle_, 'GeneralisedPrior', 'penalisation_factor',...
                value, 'f')
        end
        function value = get_penalisation_factor(self)
%***SIRF*** Returns the penalty factor in front of the prior.
            value = sirf.STIR.parameter...
                (self.handle_, 'GeneralisedPrior', 'penalisation_factor', 'f');
        end
        function set_up(self, image)
%***SIRF*** Prepares the prior for use.
            h = calllib('mstir', 'mSTIR_setupPrior',...
                self.handle_, image.handle_);
            mUtilities.check_status('Prior:set_up', h)
            mUtilities.delete(h)
        end
        function grad = get_gradient(self, image)
%***SIRF*** Returns the value of the gradient of the prior for the specified 
%         image.
            grad = sirf.STIR.ImageData();
            grad.handle_ = calllib('mstir', 'mSTIR_priorGradient', ...
                self.handle_, image.handle_);
            mUtilities.check_status('Prior:get_gradient', grad.handle_)
        end
    end
end