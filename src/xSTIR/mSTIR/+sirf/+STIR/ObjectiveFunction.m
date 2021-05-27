classdef ObjectiveFunction < handle
% Class for the iterative reconstruction algorithms objective function.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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
            name = 'ObjectiveFunction';
        end
    end
    methods
        function self = ObjectiveFunction()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_prior(self, prior)
%***SIRF*** Sets the prior (penalty term to be added to the objective function).
            sirf.Utilities.assert_validity(prior, 'Prior')
            sirf.STIR.setParameter...
                (self.handle_, 'GeneralisedObjectiveFunction', 'prior',...
                prior, 'h')
        end
        function prior = get_prior(self)
%***SIRF*** Returns the prior currently used by this objective function.
            prior = sirf.STIR.Prior();
            prior.handle_ = calllib('mstir', 'mSTIR_parameter',...
                self.handle_, 'GeneralisedObjectiveFunction', 'prior');
            sirf.Utilities.check_status...
                ('GeneralisedObjectiveFunction:get_prior', prior.handle_)
        end
        function set_num_subsets(self, num)
%***SIRF*** Sets the number of subsets of ray projections to be used 
%         for computing additive components of the gradient used by 
%         Ordered Subset algorithms for maximizing this objective function.
%         Assuming for simplicity of illustration that the ray tracing 
%         projector G is a matrix, the subsets in question correspond to
%         subsets of its rows.
            sirf.STIR.setParameter...
                (self.handle_, 'GeneralisedObjectiveFunction', ...
                'num_subsets', num, 'i')
        end
        function set_up(self, image)
%***SIRF*** Prepares this object for use.
            sirf.Utilities.assert_validity(image, 'ImageData')
            h = calllib('mstir', 'mSTIR_setupObjectiveFunction', ...
                self.handle_, image.handle_);
            sirf.Utilities.check_status('GeneralisedObjectiveFunction:set_up', h)
            sirf.Utilities.delete(h)
            %calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function v = get_value(self, image)
%***SIRF*** Returns the value of this objective function 
%         on the specified image.
            sirf.Utilities.assert_validity(image, 'ImageData')
            h = calllib('mstir', 'mSTIR_objectiveFunctionValue',...
                self.handle_, image.handle_);
            sirf.Utilities.check_status...
                ('GeneralisedObjectiveFunction:value', h)
            v = calllib('miutilities', 'mFloatDataFromHandle', h);
            sirf.Utilities.delete(h)
            %calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function g = get_subset_gradient(self, image, subset)
%***SIRF*** Returns the value of the additive component of the gradient 
%         of this objective function on the specified image corresponding 
%         to the specified subset (see method set_num_subsets()).
            sirf.Utilities.assert_validity(image, 'ImageData')
            if nargin < 3
                subset = -1;
            end
            g = sirf.STIR.ImageData();
            g.handle_ = calllib('mstir', 'mSTIR_objectiveFunctionGradient',...
                self.handle_, image.handle_, subset);
            sirf.Utilities.check_status...
                ('GeneralisedObjectiveFunction:gradient', g.handle_)
        end
        function g = get_gradient(self, image)
%***SIRF*** Returns the gradient of the objective function 
%         on the specified image.
%         image: ImageData object
            sirf.Utilities.assert_validity(image, 'ImageData')
            g = self.get_subset_gradient(image);
        end
    end
end