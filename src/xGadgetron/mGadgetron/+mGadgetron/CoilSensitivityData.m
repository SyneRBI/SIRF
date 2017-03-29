classdef CoilSensitivityData < mGadgetron.DataContainer
%     Class for a coil sensitivity maps (csm) container.
%     Each item in the container is a 4D complex array of csm values on an 
%     xyz-slice (z-dimension is normally 1).
    properties
        name_
    end
    methods (Static)
        function obj = same_object()
            obj = mGadgetron.CoilSensitivityData();
        end
    end
    methods
        function self = CoilSensitivityData()
            self.name_ = 'CoilSensitivityData';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function calculate(self, acqs)
%         Calculates coil sensitivity maps from sorted acquisitions.
%         acqs: AcquisitionData
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
%         Returns specified csm as an array of .
%         csm_num: csm (slice) number
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib...
                ('mgadgetron', 'mGT_getImageDimensions', ...
                self.handle_, csm_num - 1, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3)*dim(4);
            ptr_re = libpointer('doublePtr', zeros(n, 1));
            ptr_im = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getCoilData', ...
                self.handle_, csm_num - 1, ptr_re, ptr_im)
            re = reshape(ptr_re.Value, dim(1), dim(2), dim(3), dim(4));
            im = reshape(ptr_im.Value, dim(1), dim(2), dim(3), dim(4));
            data = complex(re, im);
        end
    end
end