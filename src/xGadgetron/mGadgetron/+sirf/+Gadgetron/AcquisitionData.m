classdef AcquisitionData < sirf.SIRF.DataContainer
% Class for MR acquisitions data.
% Each item in the container is a complex array of acquisition 
% samples for each coil.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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
    end
    methods (Static)
        function name = class_name()
            name = 'AcquisitionData';
        end
        function obj = same_object()
            obj = sirf.Gadgetron.AcquisitionData();
        end
        function set_storage_scheme(scheme)
%***SIRF*** Sets acquisition data storage scheme.
            if ~strcmp(scheme, 'memory')
                fprintf("WARNING: storage scheme '%s' not supported, ", scheme)
                fprintf('using memory storage scheme instead\n')
            end
            return
        end
        function scheme = get_storage_scheme()
%***SIRF*** Returns current acquisition storage scheme name
			scheme = 'memory';
        end
    end
    methods
        function self = AcquisitionData(filename)
%         Creates an AcquisitionData object, optionally from a raw data file
%         in HDF5 format.
            self.name_ = 'AcquisitionData';
            self.handle_ = [];
            self.sorted_ = false;
            %self.info_ = sirf.Gadgetron.AcquisitionInfo.empty(0);
            if nargin > 0
                self.handle_ = calllib('mgadgetron', ...
                    'mGT_ISMRMRDAcquisitionsFromFile', filename);
                sirf.Utilities.check_status(self.name_, self.handle_);
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function new_ad = new_acquisition_data(self, empty)
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            new_ad = sirf.Gadgetron.AcquisitionData();
            if nargin < 2
                empty = true;
            end
            if empty
                new_ad.handle_ = calllib('mgadgetron', ...
                    'mGT_createEmptyAcquisitionData', self.handle_);
            else
                new_ad.handle_ = calllib('mgadgetron', ...
                    'mGT_cloneAcquisitions', self.handle_);
            end
            sirf.Utilities.check_status(self.name_, new_ad.handle_);
        end
        function sort(self)
%***SIRF*** Sorts acquisitions with respect to (in this order):
%             - repetition
%             - phase
%             - contrast
%             - slice
%             - kspace_encode_step_2
%             - kspace_encode_step_1
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            handle = calllib('mgadgetron', 'mGT_sortAcquisitions', ...
                self.handle_);
            sirf.Utilities.check_status('AcquisitionData', handle);
            sirf.Utilities.delete(handle)
            self.sorted_ = true;
        end
        function sorted = is_sorted(self)
%***SIRF*** Returns true if acquisitions of this object are sorted
%         and false otherwise.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            sorted = (sirf.Gadgetron.parameter(self.handle_, ...
                'acquisitions', 'sorted', 'i') ~= 0);
%            sorted = self.sorted_;
        end
        function undersampled = is_undersampled(self)
%***SIRF*** Returns true if acquisitions of this object are undersampled
%         and false otherwise.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            undersampled = (sirf.Gadgetron.parameter(self.handle_, ...
                'acquisitions', 'undersampled', 'i') ~= 0);
%            sorted = self.sorted_;
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
            ap = sirf.Gadgetron.AcquisitionDataProcessor(list);
            a = ap.process(self);
        end
        function [ns, nc, na] = dimensions(self)
%***SIRF*** Returns the numbers of samples, coils and acquisitions 
%         in this AcquisitionData object.
%         Non-image-related acquisitions (noise calibration etc.) are ignored.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            ptr_i = libpointer('int32Ptr', ones(16, 1));
            calllib...
                ('mgadgetron', 'mGT_getAcquisitionDataDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
            ns = dim(1);
            nc = dim(2);
            na = prod(dim(3:end));
        end
        function a = acquisition(self, num)
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            if num < 1 || num > self.number()
                error('AcquisitionData:value_error', ...
                    'Acquisition number out of range')
            end
            a = sirf.Gadgetron.Acquisition();
            a.handle_ = calllib('mgadgetron', ...
                'mGT_acquisitionFromContainer', self.handle_, num - 1);
            sirf.Utilities.check_status('AcquisitionData', a.handle_);
        end
        function data = as_array(self)
%***SIRF*** as_array(select) returns a 3D complex array of dimensions 
%          returned by dimensions(select) containing acquisitions.
%          The meaning of select is the same as in dimensions().
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            [ns, nc, na] = self.dimensions();
            n = ns*nc*na;
            ptr_z = libpointer('singlePtr', zeros(2, n));
            calllib...
                ('mgadgetron', 'mGT_acquisitionDataAsArray', ...
                self.handle_, ptr_z, -1);
            data = reshape(ptr_z.Value(1:2:end) + 1i*ptr_z.Value(2:2:end), ...
                ns, nc, na);
        end
        function fill(self, data, select)
%***SIRF*** Changes acquisition data to that in 3D complex array argument.
            if isempty(self.handle_)
                error('AcquisitionData:empty_object', ...
                    'cannot handle empty object')
            end
            if nargin < 3
                select = 'not all';
            end
            if strcmp(select, 'all')
                all = 1;
            else
                all = 0;
            end
            z = [real(data(:))'; imag(data(:))'];
            if isa(z, 'single')
                ptr_z = libpointer('singlePtr', z);
            else
                ptr_z = libpointer('singlePtr', single(z));
            end
            h = calllib('mgadgetron', 'mGT_fillAcquisitionData', ...
                self.handle_, ptr_z, all);
            sirf.Utilities.check_status('AcquisitionData', h);
            sirf.Utilities.delete(h)
        end
    end
end
