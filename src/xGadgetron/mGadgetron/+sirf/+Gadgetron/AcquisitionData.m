classdef AcquisitionData < sirf.SIRF.DataContainer
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
%           scheme = 'file' (default):
%               all acquisition data generated from now on will be kept in
%               scratch files deleted after the user's script terminates
%           scheme = 'memory':
%               all acquisition data generated from now on will be kept in
%               RAM (avoid if data is very large)
            h = calllib...
                ('mgadgetron', 'mGT_setAcquisitionsStorageScheme', scheme);
            sirf.Utilities.check_status('AcquisitionData', h);
            sirf.Utilities.delete(h)
        end
        function scheme = get_storage_scheme()
%***SIRF*** Returns current acquisition storage scheme name
            h = calllib...
                ('mgadgetron', 'mGT_getAcquisitionsStorageScheme');
            sirf.Utilities.check_status('AcquisitionData', h);
            scheme = calllib('miutilities', 'mCharDataFromHandle', h);
            sirf.Utilities.delete(h)
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
        function sort(self)
%***SIRF*** Sorts acquisitions with respect to (in this order):
%             - repetition
%             - slice
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
            sorted = self.sorted_;
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
            a = sirf.Gadgetron.Acquisition();
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
            [ns, nc, na] = self.dimensions(select);
            %na = self.number();
            if strcmp(select, 'all')
                all = 1;
            else
                all = 0;
            end
            n = ns*nc*na;
            ptr_z = libpointer('singlePtr', zeros(2, n));
            calllib...
                ('mgadgetron', 'mGT_acquisitionsDataAsArray', ...
                self.handle_, ptr_z, all);
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
                select = 'all';
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
            h = calllib('mgadgetron', 'mGT_fillAcquisitionsData', ...
                self.handle_, ptr_z, all);
            sirf.Utilities.check_status('AcquisitionData', h);
            sirf.Utilities.delete(h)
        end
    end
end
