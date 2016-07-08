classdef PoissonLogLh_LinModMean_AcqModData < stir.PoissonLogLh_LinModMean
    methods
        function self = PoissonLogLh_LinModMean_AcqModData(obj_fun)
            self.name =...
                'PoissonLogLikelihoodWithLinearModelForMeanAndProjData';
            if nargin < 1
                self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            else
                self.handle = calllib...
                    ('mutilities', 'mCopyOfObjectHandle', obj_fun.handle);
%                    ('mstir', 'mSTIR_copyOfObject', obj_fun.handle);
            end
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
%            calllib('mstir', 'mSTIR_deleteObject', self.handle)
            self.handle = [];
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
        function set_acquisition_model(self, am)
            stir.setParameter(self.handle, self.name,...
                'projector_pair_type', am.handle, 'h')
        end
        function am = get_acquisition_model(self)
            am = stir.AcquisitionModelUsingMatrix();
            am.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'projector_pair_type');
            stir.checkExecutionStatus...
                ([self.name ':get_projector_pair'], am.handle)
        end
        function set_acquisition_data(self, am)
            stir.setParameter(self.handle, self.name,...
                'proj_data_sptr', am.handle, 'h')
        end
    end
end