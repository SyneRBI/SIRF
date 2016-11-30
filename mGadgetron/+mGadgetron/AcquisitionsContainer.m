classdef AcquisitionsContainer < mGadgetron.DataContainer
    properties
        sorted_
    end
    methods
        function self = AcquisitionsContainer()
            self.handle_ = [];
            self.sorted_ = false;
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function sort(self)
            handle = calllib('mgadgetron', 'mGT_orderAcquisitions', ...
                self.handle_);
            mGadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            self.sorted_ = true;
        end
        function sorted = is_sorted(self)
            sorted = self.sorted_;
        end
        function a = process(self, list)
            ap = mGadgetron.AcquisitionsProcessor(list);
            a = ap.process(self);
        end
        function [ns, ny, nc] = slice_dimensions(self)
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            ns = dim(1);
            ny = dim(2);
            nc = dim(3);
        end
        function data = slice_as_array(self, num)
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
            ptr_re = libpointer('doublePtr', zeros(n, 1));
            ptr_im = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsData', ...
                self.handle_, num - 1, ptr_re, ptr_im);
            re = reshape(ptr_re.Value, dim(1), dim(2), dim(3));
            im = reshape(ptr_im.Value, dim(1), dim(2), dim(3));
            data = re + 1i*im;
        end
    end
end