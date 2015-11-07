classdef PoissonLogLikelihoodWithLinearModelForMeanAndProjData...
        < stir.PoissonLogLikelihoodWithLinearModelForMean
    properties
        projectors
    end
    methods
        function self =...
                PoissonLogLikelihoodWithLinearModelForMeanAndProjData(obj_fun)
            self.name =...
                'PoissonLogLikelihoodWithLinearModelForMeanAndProjData';
            if nargin < 1
                self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
                self.owns_handle = true;
            else
                self.handle = calllib('mstir', 'mRefDataHandle', obj_fun.handle);
                self.owns_handle = false;
            end
            self.projectors = [];
        end
        function delete(self)
            if self.owns_handle
                calllib('mstir', 'mSTIR_deleteObject', self.handle,...
                    'ObjectiveFunction')
                self.handle = [];
            end
        end
        function set_input_filename(self, filename)
            stir.setParameter(self.handle, self.name,...
                'input_filename', filename, 'c')
        end
        function set_zero_seg0_end_planes(self, flag)
            if flag
                str = 'true';
            else
                str = 'false';
            end
            stir.setParameter(self.handle, self.name,...
                'zero_seg0_end_planes', str, 'c') 
        end
        function set_max_segment_num_to_process(self, n)
            stir.setParameter(self.handle, self.name, ...
                'max_segment_num_to_process', n, 'i') 
        end
        function set_projector_pair(self, pp)
            stir.setParameter(self.handle, self.name,...
                'projector_pair_type', pp.handle, 'h')
            self.projectors = pp;
        end
        function projectors = get_projector_pair(self)
            projectors = self.projectors;
            if isempty(projectors)
                error([self.name ':no_projectors'],...
                    [self.name ': no projectors set'])
            end
        end
    end
end