classdef IterativeReconstruction < stir.Reconstruction
    properties (Constant)
        IR = 'IterativeReconstruction';
    end
    methods
        function self = IterativeReconstruction()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
%                calllib('mstir', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_num_subsets(self, n)
            stir.setParameter(self.handle, self.IR, 'num_subsets', n, 'i')
        end
        function n = get_num_subsets(self)
            n = stir.parameter(self.handle, self.IR, 'num_subsets', 'i');
        end
        function set_start_subset_num(self, n)
            stir.setParameter(self.handle,...
                self.IR, 'start_subset_num', n, 'i')
        end
        function n = get_start_subset_num(self)
            n = stir.parameter(self.handle,...
                self.IR, 'start_subset_num', 'i');
        end
        function set_num_subiterations(self, n)
            stir.setParameter(self.handle,...
                self.IR, 'num_subiterations', n, 'i')
        end
        function n = get_num_subiterations(self)
            n = stir.parameter(self.handle,...
                self.IR, 'num_subiterations', 'i');
        end
        function set_start_subiteration_num(self, n)
            stir.setParameter(self.handle,...
                self.IR, 'start_subiteration_num', n, 'i')
        end
        function n = get_start_subiteration_num(self)
            n = stir.parameter(self.handle,...
                self.IR, 'start_subiteration_num', 'i');
        end
        function set_subiteration_num(self, iter)
            stir.setParameter(self.handle,...
                self.IR, 'subiteration_num', iter, 'i')
        end
        function iter = get_subiteration_num(self)
            iter = stir.parameter(self.handle,...
                self.IR, 'subiteration_num', 'i');
        end
        function set_save_interval(self, n)
            stir.setParameter(self.handle,...
                self.IR, 'save_interval', n, 'i')
        end
        function set_inter_iteration_filter_interval(self, n)
            stir.setParameter(self.handle, self.IR,...
                'inter_iteration_filter_interval', n, 'i')
        end
        function set_objective_function(self, obj_fun)
            stir.setParameter(self.handle, self.IR,...
                'objective_function', obj_fun.handle, 'h')
        end
        function set_inter_iteration_filter(self, filter)
            stir.setParameter(self.handle, self.IR,...
                'inter_iteration_filter_type', filter.handle, 'h')
        end
        function filter = get_inter_iteration_filter(self)
            filter = stir.DataProcessor();
            filter.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.IR, 'inter_iteration_filter_type');
            stir.checkExecutionStatus...
                ([self.IR ':get_inter_iteration_filter'], filter.handle)
        end
        function set_up(self, image)
            h = calllib('mstir', 'mSTIR_setupReconstruction',...
                self.handle, image.handle);
            stir.checkExecutionStatus([self.IR ':set_up'], h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function update(self, image)
            h = calllib('mstir', 'mSTIR_updateReconstruction',...
                self.handle, image.handle);
            stir.checkExecutionStatus([self.IR ':update'], h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function reconstruct(self, image)
            h = calllib('mstir', 'mSTIR_runReconstruction',...
                self.handle, image.handle);
            stir.checkExecutionStatus([self.IR ':reconstruct'], h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
    end
end