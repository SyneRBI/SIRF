classdef IterativeReconstructor < mStir.Reconstructor
% Class for a generic iterative PET reconstructor.

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

    properties (Constant)
        IR = 'IterativeReconstruction';
    end
    properties
        input
        image
        subset
    end
    methods
        function self = IterativeReconstructor()
            self.handle = [];
            self.input = [];
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
            mStir.setParameter(self.handle, self.IR, 'num_subsets', n, 'i')
        end
        function n = get_num_subsets(self)
            n = mStir.parameter(self.handle, self.IR, 'num_subsets', 'i');
        end
        function set_start_subset_num(self, n)
            mStir.setParameter(self.handle,...
                self.IR, 'start_subset_num', n, 'i')
        end
        function n = get_start_subset_num(self)
            n = mStir.parameter(self.handle,...
                self.IR, 'start_subset_num', 'i');
        end
        function set_num_subiterations(self, n)
            mStir.setParameter(self.handle,...
                self.IR, 'num_subiterations', n, 'i')
        end
        function n = get_num_subiterations(self)
            n = mStir.parameter(self.handle,...
                self.IR, 'num_subiterations', 'i');
        end
        function set_start_subiteration_num(self, n)
            mStir.setParameter(self.handle,...
                self.IR, 'start_subiteration_num', n, 'i')
        end
        function n = get_start_subiteration_num(self)
            n = mStir.parameter(self.handle,...
                self.IR, 'start_subiteration_num', 'i');
        end
        function set_subiteration_num(self, iter)
            mStir.setParameter(self.handle,...
                self.IR, 'subiteration_num', iter, 'i')
        end
        function iter = get_subiteration_num(self)
            iter = mStir.parameter(self.handle,...
                self.IR, 'subiteration_num', 'i');
        end
        function set_save_interval(self, n)
            mStir.setParameter(self.handle,...
                self.IR, 'save_interval', n, 'i')
        end
        function set_inter_iteration_filter_interval(self, n)
            mStir.setParameter(self.handle, self.IR,...
                'inter_iteration_filter_interval', n, 'i')
        end
        function set_objective_function(self, obj_fun)
            mStir.setParameter(self.handle, self.IR,...
                'objective_function', obj_fun, 'h')
        end
        function set_inter_iteration_filter(self, filter)
            mStir.setParameter(self.handle, self.IR,...
                'inter_iteration_filter_type', filter, 'h')
        end
        function filter = get_inter_iteration_filter(self)
            filter = mStir.DataProcessor();
            filter.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.IR, 'inter_iteration_filter_type');
            mUtil.checkExecutionStatus...
                ([self.IR ':get_inter_iteration_filter'], filter.handle)
        end
        function set_up(self, image)
            h = calllib('mstir', 'mSTIR_setupReconstruction',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus([self.IR ':set_up'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function set_input(self, input)
            self.input = input;
        end
        function set_current_subset(self, subset)
            self.subset = subset;
        end
        function set_current_estimate(self, image)
            self.image = image;
        end
        function image = get_current_estimate(self)
            image = self.image;
        end
        function update_current_estimate(self)
            if isempty(self.image)
                error([self.IR ':update_current_image'], ...
                    'current estimate not set')
            end
            h = calllib('mstir', 'mSTIR_updateReconstruction',...
                self.handle, self.image.handle);
            mUtil.checkExecutionStatus([self.IR ':update_current_image'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function process(self)
            if isempty(self.image)
                error([self.IR ':process'], 'current estimate not set')
            end
            h = calllib('mstir', 'mSTIR_runReconstruction',...
                self.handle, self.image.handle);
            mUtil.checkExecutionStatus([self.IR ':process'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function update(self, image)
            h = calllib('mstir', 'mSTIR_updateReconstruction',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus([self.IR ':update'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function reconstruct(self, image)
            h = calllib('mstir', 'mSTIR_runReconstruction',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus([self.IR ':reconstruct'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
    end
end