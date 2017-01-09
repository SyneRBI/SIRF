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
        function [ns, nc, na] = dimensions(self, select)
            ptr_i = libpointer('int32Ptr', ones(16, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            all = true;
            if nargin > 1
                all = strcmp(select, 'all');
            end
            ns = dim(1);
            nc = dim(2);
            if all
                na = self.number();
            else
                na = prod(dim(3:end));
            end
        end
        function data = as_array(self, select)
            if nargin < 2
                select = 'all';
            end
            [ns, nc, ma] = self.dimensions(select);
            na = self.number();
            if strcmp(select, 'all')
                n = na;
                ma = na;
            else
                n = na + 1;
            end
            m = ns*nc*ma;
            ptr_re = libpointer('doublePtr', zeros(m, 1));
            ptr_im = libpointer('doublePtr', zeros(m, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsData', ...
                self.handle_, n, ptr_re, ptr_im);
            re = reshape(ptr_re.Value, ns, nc, ma);
            im = reshape(ptr_im.Value, ns, nc, ma);
            data = re + 1i*im;
        end
        function fill(self, data)
            [ns, nc, na] = size(data);
            re = real(data);
            im = imag(data);
            ptr_re = libpointer('doublePtr', re);
            ptr_im = libpointer('doublePtr', im);
            h = calllib('mgadgetron', 'mGT_setAcquisitionsData', ...
                self.handle_, na, nc, ns, ptr_re, ptr_im);
            mGadgetron.checkExecutionStatus('AcquisitionsContainer', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function [ns, ny, nc] = slice_dimensions(self)
            ptr_i = libpointer('int32Ptr', ones(16, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            ns = dim(1);
            nc = dim(2);
            ny = dim(3);
        end
        function data = slice_as_array(self, num)
            ptr_i = libpointer('int32Ptr', ones(16, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            ns = dim(1);
            nc = dim(2);
            ny = dim(3);
            n = ns*nc*ny;
            %n = dim(1)*dim(2)*dim(3);
            ptr_re = libpointer('doublePtr', zeros(n, 1));
            ptr_im = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsData', ...
                self.handle_, num - 1, ptr_re, ptr_im);
            re = reshape(ptr_re.Value, ns, ny, nc);
            im = reshape(ptr_im.Value, ns, ny, nc);
%             re = reshape(ptr_re.Value, dim(1), dim(2), dim(3));
%             im = reshape(ptr_im.Value, dim(1), dim(2), dim(3));
            data = re + 1i*im;
        end
    end
end