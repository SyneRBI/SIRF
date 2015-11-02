classdef IterativeReconstruction < stir.Reconstruction
    properties
        obj_fun
    end
    methods
        function self = IterativeReconstruction()
            self.obj_fun = [];
        end
        function set_num_subsets(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'num_subsets', n, 'i')
        end
        function n = get_num_subsets(self)
            n = stir.parameter(self.handle,...
                'IterativeReconstruction', 'num_subsets', 'i');
        end
        function set_start_subset_num(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'start_subset_num', n, 'i')
        end
        function n = get_start_subset_num(self)
            n = stir.parameter(self.handle,...
                'IterativeReconstruction', 'start_subset_num', 'i');
        end
        function set_num_subiterations(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'num_subiterations', n, 'i')
        end
        function n = get_num_subiterations(self)
            n = stir.parameter(self.handle,...
                'IterativeReconstruction', 'num_subiterations', 'i');
        end
        function set_start_subiteration_num(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'start_subiteration_num', n, 'i')
        end
        function n = get_start_subiteration_num(self)
            n = stir.parameter(self.handle,...
                'IterativeReconstruction', 'start_subiteration_num', 'i');
        end
        function set_subiteration_num(self, iter)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'subiteration_num', iter, 'i')
        end
        function iter = get_subiteration_num(self)
            iter = stir.parameter(self.handle,...
                'IterativeReconstruction', 'subiteration_num', 'i');
        end
        function set_save_interval(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction', 'save_interval', n, 'i')
        end
        function set_inter_iteration_filter_interval(self, n)
            stir.setParameter(self.handle,...
                'IterativeReconstruction',...
                'inter_iteration_filter_interval', n, 'i')
        end
        function set_objective_function(self, obj_fun)
            stir.setParameter(self.handle,...
                'IterativeReconstruction',...
                'objective_function', obj_fun.handle, 'h')
            self.obj_fun = obj_fun;
        end
        function obj_fun = get_objective_function(self)
            obj_fun = self.obj_fun;
            if isempty(obj_fun)
                error('IterativeReconstruction:no_obj_fun_set',...
                    'IterativeReconstruction: no objective function set')
            end
        end
        function set_inter_iteration_filter(self, filter)
            stir.setParameter(self.handle,...
                'IterativeReconstruction',...
                'inter_iteration_filter_type', filter.handle, 'h')
        end
        function set_up(self, image)
            h = calllib('mstir', 'mSTIR_setupReconstruction',...
                self.handle, image.handle);
            stir.checkExecutionStatus('set_up', h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function update(self, image)
            h = calllib('mstir', 'mSTIR_update', self.handle, image.handle);
            stir.checkExecutionStatus('update', h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function reconstruct(self, image)
            h = calllib('mstir', 'mSTIR_reconstruct',...
                self.handle, image.handle);
            stir.checkExecutionStatus('reconstruct', h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
    end
end