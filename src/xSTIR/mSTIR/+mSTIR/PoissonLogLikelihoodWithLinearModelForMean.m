classdef PoissonLogLikelihoodWithLinearModelForMean < mSTIR.ObjectiveFunction
% ADVANCED USERS ONLY.
% Class for STIR PoissonLogLikelihoodWithLinearModelForMean object, see
% http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMean.html

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

    methods
        function self = PoissonLogLikelihoodWithLinearModelForMean()
%         Creates new empty object.
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_sensitivity_filename(self, name)
%***SIRF*** Specifies the file with the sensitivity data to be used.
            mSTIR.setParameter...
                (self.handle_, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'sensitivity_filename', name, 'c')
        end
        function set_use_subset_sensitivities(self, value)
%***SIRF*** Specifies whether the subset sensitivity data is to be used.
            if value
                str = 'true';
            else
                str = 'false';
            end
            mSTIR.setParameter...
                (self.handle_, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'use_subset_sensitivities', str, 'c')
        end
        function set_recompute_sensitivity(self, value)
%***SIRF*** Specifies whether the subset sensitivity data must be re-computed.
            if value
                str = 'true';
            else
                str = 'false';
            end
            mSTIR.setParameter...
                (self.handle_, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'recompute_sensitivity', str, 'c')
        end
        function sens = get_subset_sensitivity(self, subset)
%***SIRF*** Returns the specified subset sensitivity data as ImageData.
            sens = mSTIR.ImageData();
            sens.handle_ = calllib('mstir', 'mSTIR_subsetSensitivity',...
                self.handle_, subset);
            mUtilities.check_status...
                ('PoissonLinModMean:get_subset_sensitivity',...
                sens.handle_)
        end
        function bar = get_backprojection_of_acquisition_ratio...
                (self, image, subset)
%***SIRF*** Returns the backprojection of the ratio of measured to estimated
%         acquisition data for the specified image and subset.
            mUtilities.assert_validity(image, 'ImageData')
            bar = mSTIR.ImageData();
            bar.handle_ = calllib...
                ('mstir', 'mSTIR_objectiveFunctionGradientNotDivided',...
                self.handle_, image.handle_, subset);
            mUtilities.check_status...
                ('PoissonLinModMean:get_backprojection_of_acquisition_ratio',...
                bar.handle_)
        end
    end
end
