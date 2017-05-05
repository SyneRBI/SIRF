classdef AcquisitionData < mGadgetron.DataContainer
% Class for MR acquisitions data.
% Each item in the container is a complex array of acquisition 
% samples for each coil.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

    properties
        % Class name
        name_
        % Are acquisitions sorted?
        sorted_
        % array of AcquisitionInfo objects
        info_
    end
    methods (Static)
        function obj = same_object()
            obj = mGadgetron.AcquisitionData();
        end
    end
    methods
        function self = AcquisitionData(filename)
%         Creates an AcquisitionData object, optionally from a raw data file
%         in HDF5 format.
            self.name_ = 'AcquisitionData';
            self.handle_ = [];
            self.sorted_ = false;
            self.info_ = mGadgetron.AcquisitionInfo.empty(0);
            if nargin > 0
                self.handle_ = calllib('mgadgetron', ...
                    'mGT_ISMRMRDAcquisitionsFromFile', filename);
                mUtil.checkExecutionStatus(self.name_, self.handle_);
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function sort(self)
%***SIRF*** Sorts acquisitions with respect to (in this order):
%             - repetition
%             - slice
%             - kspace_encode_step_1
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            handle = calllib('mgadgetron', 'mGT_orderAcquisitions', ...
                self.handle_);
            mUtil.checkExecutionStatus('AcquisitionData', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            self.sorted_ = true;
        end
        function sorted = is_sorted(self)
%***SIRF*** Returns true if acquisitions of this object are sorted
%         and false otherwise.
            sorted = self.sorted_;
        end
        function set_info(self)
%***SIRF*** Fills the acquisition info array (a property of this class)
%         with information on the data acquisition (e.g. type of data 
%         (noise, image,..), k-space location, slice number,...)           
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
        function info = get_info(self, info_type)
%***SIRF*** Returns an array of acquisitions info of the specified type.
%         Currently implemented types are: 
%         - flags
%         - encode_step_1
%         - slice
%         - repetition
            if isempty(self.info_)
                self.set_info()
            end
            na = self.number();
            info = zeros(na, 1);
            if strcmp(info_type, 'flags')
                for a = 1 : na
                    info(a) = self.info_(a).flags_;
                end
            elseif strcmp(info_type, 'encode_step_1')
                for a = 1 : na
                    info(a) = self.info_(a).encode_step_1_;
                end
            elseif strcmp(info_type, 'slice')
                for a = 1 : na
                    info(a) = self.info_(a).slice_;
                end
            elseif strcmp(info_type, 'repetition')
                for a = 1 : na
                    info(a) = self.info_(a).repetition_;
                end
            end
        end
        function a = process(self, list)
%***SIRF*** Returns acquisitions processed by a chain of gadgets.
%         The argument is a cell array of gadget definitions
%         [{gadget1_definition}, {gadget2_definition}, ...],
%         where gadget definitions are strings of the form
%         'label:name(property1=value1,property2=value2,...)',
%         where the only mandatory field is name, the Gadgetron 
%         name of the gadget. An optional expression in round 
%         brackets can be used to assign values to gadget properties,
%         and an optional label can be used to change the labelled
%         gadget properties after the chain has been defined.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            ap = mGadgetron.AcquisitionDataProcessor(list);
            a = ap.process(self);
        end
        function a = clone(self)
%***SIRF*** Returns a copy of self.
            ap = mGadgetron.AcquisitionDataProcessor();
            a = ap.process(self);
        end
        function [ns, nc, na] = dimensions(self, select)
%***SIRF*** Returns the numbers of samples, coils and acquisitions 
%         in this AcquisitionData object.
%         If the argument is supplied that is not 'all', then non-image 
%         related acquisitions (noise calibration etc.) are ignored.
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
%***SIRF*** as_array(select) returns an array with this object's data 
%         (a 3D complex array).
%         The dimensions are those returned by dimensions(select).
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
%***SIRF*** Changes acquisition data to that in 3D complex array argument.
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
            mUtil.checkExecutionStatus('AcquisitionData', h);
            calllib('mutilities', 'mDeleteDataHandle', self.handle_)
            self.handle_ = h;
        end
    end
end
