classdef PoissonLogLikelihoodWithLinearModelForMeanAndProjData...
        < stir.PoissonLogLikelihoodWithLinearModelForMean
    properties
        projectors
    end
    methods
        function self = PoissonLogLikelihoodWithLinearModelForMeanAndProjData()
            %fprintf('derived class constructor called\n')
            self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData';
            calllib('mstir', 'mDeleteDataHandle', self.handle)
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            self.projectors = [];
        end
        function delete(self)
            %fprintf('derived class destructor called\n')
            calllib('mstir', 'mSTIR_deleteObject', self.handle, self.name)
            self.handle = calllib('mstir', 'mNewDataHandle');
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