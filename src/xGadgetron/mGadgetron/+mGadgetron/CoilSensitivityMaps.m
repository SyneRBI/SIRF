classdef CoilSensitivityMaps < mGadgetron.DataContainer
    properties
        name_
    end
    methods
        function self = CoilSensitivityMaps()
            self.name_ = 'CoilSensitivityMaps';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function calculate(self, acqs)
            if ~acqs.is_sorted()
                fprintf('WARNING: acquisitions may be in a wrong order\n')
            end
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = calllib('mgadgetron', 'mGT_CoilSensitivities', '');
            mUtil.checkExecutionStatus(self.name_, self.handle_);
            handle = calllib('mgadgetron', 'mGT_computeCoilSensitivities', ...
                self.handle_, acqs.handle_);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function data = csm_as_array(self, csm_num)
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib...
                ('mgadgetron', 'mGT_getImageDimensions', ...
                self.handle_, csm_num - 1, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3)*dim(4);
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getCMSDataAbs', ...
                self.handle_, csm_num - 1, ptr_v)
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3), dim(4));
        end
    end
end