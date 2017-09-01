classdef IterativeReconstructor < mSTIR.Reconstructor
% Class for a generic iterative PET reconstructor.
% Iterative reconstructors of this class maximize the objective function
% that is formed by components corresponding to subsets of the acquisition
% data.

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

    properties (Constant)
        IR = 'IterativeReconstruction';
    end
    properties
        subset
    end
    methods
        function self = IterativeReconstructor()
            self.handle = [];
            self.image = [];
            self.subset = 0;
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_num_subsets(self, n)
%***SIRF*** Sets the number of components of acquisition data.
            mSTIR.setParameter(self.handle, self.IR, 'num_subsets', n, 'i')
        end
        function n = get_num_subsets(self)
%***SIRF*** Returns the number of components of acquisition data.
            n = mSTIR.parameter(self.handle, self.IR, 'num_subsets', 'i');
        end
        function set_start_subset_num(self, n)
%***SIRF*** Sets the first component of acquisition data to work with.
            mSTIR.setParameter(self.handle,...
                self.IR, 'start_subset_num', n, 'i')
        end
        function n = get_start_subset_num(self)
%***SIRF*** Returns the first component of acquisition data to work with.
            n = mSTIR.parameter(self.handle,...
                self.IR, 'start_subset_num', 'i');
        end
        function set_num_subiterations(self, n)
%***SIRF*** Sets the number of subset iterations.
            mSTIR.setParameter(self.handle,...
                self.IR, 'num_subiterations', n, 'i')
        end
        function n = get_num_subiterations(self)
%***SIRF*** Returns the number of subset iterations.
            n = mSTIR.parameter(self.handle,...
                self.IR, 'num_subiterations', 'i');
        end
        function set_start_subiteration_num(self, n)
%***SIRF*** Sets the subset iteration to start with.
            mSTIR.setParameter(self.handle,...
                self.IR, 'start_subiteration_num', n, 'i')
        end
        function n = get_start_subiteration_num(self)
%***SIRF*** Returns the subset iteration to start with.
            n = mSTIR.parameter(self.handle,...
                self.IR, 'start_subiteration_num', 'i');
        end
        function set_subiteration_num(self, iter)
%***SIRF*** Sets the current subset iteration.
            mSTIR.setParameter(self.handle,...
                self.IR, 'subiteration_num', iter, 'i')
        end
        function iter = get_subiteration_num(self)
%***SIRF*** Returns the current subset iteration.
            iter = mSTIR.parameter(self.handle,...
                self.IR, 'subiteration_num', 'i');
        end
        function set_save_interval(self, n)
%***SIRF*** Sets the saving interval size. 
%         Size 1 means the iterated approximate image must be saved to 
%         a file every iteration, size 2 every other iteration etc.
            mSTIR.setParameter(self.handle,...
                self.IR, 'save_interval', n, 'i')
        end
%         function set_inter_iteration_filter_interval(self, n)
% %***SIRF*** Sets the filtering interval size: size 1 means the iterated 
% %         approximate image must be filtered every iteration, size 2 every 
% %         other iteration etc.
%             mSTIR.setParameter(self.handle, self.IR,...
%                 'inter_iteration_filter_interval', n, 'i')
%         end
        function set_objective_function(self, obj_fun)
%***SIRF*** Sets the objective function to be maximized.
            mSTIR.setParameter(self.handle, self.IR,...
                'objective_function', obj_fun, 'h')
        end
%         function set_inter_iteration_filter(self, filter)
% %***SIRF*** Sets the filter to be applied at intervals specified by current
% %         setting for inter-iteration filter interval.
%             mSTIR.setParameter(self.handle, self.IR,...
%                 'inter_iteration_filter_type', filter, 'h')
%         end
%         function filter = get_inter_iteration_filter(self)
% %***SIRF*** Returns the inter-iteration filter currently in use.
%             filter = mSTIR.DataProcessor();
%             filter.handle = calllib('mstir', 'mSTIR_parameter',...
%                 self.handle, self.IR, 'inter_iteration_filter_type');
%             mUtilities.check_status...
%                 ([self.IR ':get_inter_iteration_filter'], filter.handle)
%         end
        function set_up(self, image)
%***SIRF*** Prepares the reconstructor for use.
%         The argumant is an ImageData object used as a template for the
%         reconstructed image.
            h = calllib('mstir', 'mSTIR_setupReconstruction',...
                self.handle, image.handle);
            mUtilities.check_status([self.IR ':set_up'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function set_current_subset(self, subset)
%***SIRF*** Sets the current subset to work with.
            self.subset = subset;
        end
        function set_current_estimate(self, image)
%***SIRF*** Sets the current image estimate.
            self.image = image;
        end
        function image = get_current_estimate(self)
%***SIRF*** Returns the current image estimate.
            image = self.image;
        end
        function update_current_estimate(self)
%***SIRF*** Updates the current image estimate.
%         This uses data from the current subset to update the current 
%         estimate, i.e. it performs one sub-iteration.
            if isempty(self.image)
                error([self.IR ':update_current_image'], ...
                    'current estimate not set')
            end
            h = calllib('mstir', 'mSTIR_updateReconstruction',...
                self.handle, self.image.handle);
            mUtilities.check_status([self.IR ':update_current_image'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
%         function process(self)
% %***SIRF*** Reconstruct the image 
% %         by applying currently set range of
% %         iterations to the current image estimate.
%             if isempty(self.image)
%                 error([self.IR ':process'], 'current estimate not set')
%             end
%             h = calllib('mstir', 'mSTIR_runReconstruction',...
%                 self.handle, self.image.handle);
%             mUtilities.check_status([self.IR ':process'], h)
%             calllib('mutilities', 'mDeleteDataHandle', h)
%         end
        function update(self, image)
%***SIRF*** Updates the image estimate specified by the argument 
%         by performing one iteration on the current subspace.
            h = calllib('mstir', 'mSTIR_updateReconstruction',...
                self.handle, image.handle);
            mUtilities.check_status([self.IR ':update'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
%         function reconstruct(self, image)
% %***SIRF*** Reconstruct the image 
% %         by applying currently set range of
% %         iterations to the image estimate specified by the argument.
%             h = calllib('mstir', 'mSTIR_runReconstruction',...
%                 self.handle, image.handle);
%             mUtilities.check_status([self.IR ':reconstruct'], h)
%             calllib('mutilities', 'mDeleteDataHandle', h)
%         end
    end
end