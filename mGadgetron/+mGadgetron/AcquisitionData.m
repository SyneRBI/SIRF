classdef AcquisitionData < mGadgetron.DataContainer %mGadgetron.AcquisitionsContainer
    % Class for MR acquisition data
    properties
        % Class name
        name_
        % Are acquisitions sorted?
        sorted_
        % array of AcquisitionInfo objects
        info_
    end
    methods
        function self = AcquisitionData(filename)
            self.name_ = 'AcquisitionData';
            self.handle_ = [];
            self.sorted_ = false;
            self.info_ = mGadgetron.AcquisitionInfo.empty(0);
            if nargin > 0
                self.handle_ = calllib('mgadgetron', ...
                    'mGT_ISMRMRDAcquisitionsFromFile', filename);
                mGadgetron.checkExecutionStatus(self.name_, self.handle_);
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function sort(self)
            % Sorts acquisitions.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            handle = calllib('mgadgetron', 'mGT_orderAcquisitions', ...
                self.handle_);
            mGadgetron.checkExecutionStatus('AcquisitionData', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            self.sorted_ = true;
        end
        function sorted = is_sorted(self)
            sorted = self.sorted_;
        end
        function set_info(self)
            na = self.number();
            self.info_(na) = mGadgetron.AcquisitionInfo;
            for ia = 1 : na
                acq = self.acquisition(ia);
                info = mGadgetron.AcquisitionInfo();
                info.flags_ = acq.flags();
                info.encode_step_1_ = acq.idx_kspace_encode_step_1();
                info.slice_ = acq.idx_slice();
                info.repetition_ = acq.idx_repetition();
                self.info_(ia) = info;
            end
        end
        function info = get_info(self, par)
            if isempty(self.info_)
                self.set_info()
            end
            na = self.number();
            info = zeros(na);
            if strcmp(par, 'flags')
                for a = 1 : na
                    info(a) = self.info_(a).flags_;
                end
            elseif strcmp(par, 'encode_step_1')
                for a = 1 : na
                    info(a) = self.info_(a).encode_step_1_;
                end
            elseif strcmp(par, 'slice')
                for a = 1 : na
                    info(a) = self.info_(a).slice_;
                end
            elseif strcmp(par, 'repetition')
                for a = 1 : na
                    info(a) = self.info_(a).repetition_;
                end
            end
        end
        function a = process(self, list)
            % Returns acquisitions processed by a chain of gadgets.
            % The argument is a cell array of gadget definitions
            % [{gadget1_definition}, {gadget2_definition}, ...],
            % where gadget definitions are strings of the form
            % 'label:name(property1=value1,property2=value2,...)',
            % where the only mandatory field is name, the Gadgetron 
            % name of the gadget. An optional expression in round 
            % brackets can be used to assign values to gadget properties,
            % and an optional label can be used to change the labelled
            % gadget properties after the chain has been defined.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            ap = mGadgetron.AcquisitionsProcessor(list);
            a = ap.process(self);
        end
        function [ns, nc, na] = dimensions(self, select)
            % Finds number of samples, coils and acquisitions.
            % If the argument is supplied that is not 'all', then
            % non-image related acquisitions are ignored.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
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
        function a = acquisition(self, num)
            a = mGadgetron.Acquisition();
            a.handle_ = calllib('mgadgetron', ...
                'mGT_acquisitionFromContainer', self.handle_, num - 1);
        end
        function data = as_array(self, select)
            % Returns 3D complex array containing acquisition data.
            % The dimensions are those returned by dimensions().
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
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
            % Changes acquisition data to that provided by 3D array argument.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            [ns, nc, na] = size(data);
            re = real(data);
            im = imag(data);
            ptr_re = libpointer('doublePtr', re);
            ptr_im = libpointer('doublePtr', im);
            h = calllib('mgadgetron', 'mGT_setAcquisitionsData', ...
                self.handle_, na, nc, ns, ptr_re, ptr_im);
            mGadgetron.checkExecutionStatus('AcquisitionData', h);
            calllib('mutilities', 'mDeleteDataHandle', self.handle_)
            self.handle_ = h;
        end
        function [ns, ny, nc] = slice_dimensions(self)
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
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
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            ptr_i = libpointer('int32Ptr', ones(16, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            ns = dim(1);
            nc = dim(2);
            ny = dim(3);
            n = ns*nc*ny;
            ptr_re = libpointer('doublePtr', zeros(n, 1));
            ptr_im = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionsData', ...
                self.handle_, num - 1, ptr_re, ptr_im);
            re = reshape(ptr_re.Value, ns, ny, nc);
            im = reshape(ptr_im.Value, ns, ny, nc);
            data = re + 1i*im;
        end
    end
end