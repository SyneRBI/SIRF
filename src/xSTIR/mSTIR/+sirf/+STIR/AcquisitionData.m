classdef AcquisitionData < sirf.SIRF.DataContainer
% Class for PET acquisition data objects.

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
        name
        read_only
    end
    methods (Static)
        function name = class_name()
            name = 'AcquisitionData';
        end
        function obj = same_object()
            obj = sirf.STIR.AcquisitionData();
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
                ('mstir', 'mSTIR_setAcquisitionDataStorageScheme', scheme);
            sirf.Utilities.check_status('AcquisitionData', h);
            sirf.Utilities.delete(h)
        end
        function scheme = get_storage_scheme()
%***SIRF*** Returns current acquisition storage scheme name
            h = calllib...
                ('mstir', 'mSTIR_getAcquisitionDataStorageScheme');
            sirf.Utilities.check_status('AcquisitionData', h);
            scheme = calllib('miutilities', 'mCharDataFromHandle', h);
            sirf.Utilities.delete(h)
        end
    end
    methods
        function self = AcquisitionData...
                (arg, span, max_ring_diff, view_mash_factor)
%***SIRF*** AcquisitionData(arg) creates a new AcquisitionData object 
%           from a file or scanner or another AcquisitionData object;
%           arg:  file or scanner name or AcquisitionData object.
%           if a scanner name is used, additional arguments can be
%           given to specify the data size, e.g.:
%                acq=AcquisitionData('Siemens_mMR',span,max_ring_diff,view_mash_factor);
%           Defaults are: 
%                span=1, max_ring_diff=-1 (i.e. all), view_mash_factor=1
            self.handle_ = [];
            self.name = 'AcquisitionData';
            self.read_only = false;
            if nargin < 1
                return
            elseif ischar(arg)
                i = strfind(arg, '.');
                if isempty(i)
                    if nargin < 4
                        view_mash_factor = 1;
                    end
                    if nargin < 3
                        max_ring_diff = -1;
                    end
                    if nargin < 2
                        span = 1;
                    end
                    self.handle_ = calllib...
                        ('mstir', 'mSTIR_acquisitionDataFromScannerInfo',...
                        arg, span, max_ring_diff, view_mash_factor);
                    status = calllib('miutilities', 'mExecutionStatus', ...
                        self.handle_);
                    if status ~= 0
                        msg = calllib('miutilities', 'mExecutionError', ...
                            self.handle_);
                        if strcmp(msg, 'Unknown scanner')
                            err_msg_fmt = ['Unknown scanner %s or ' ...
                                'missing raw data file extension'];
                            error('AcquisitionData:wrong_data_source', ...
                                err_msg_fmt, arg)
                        end
                    end
                else
                    self.handle_ = calllib('mstir', 'mSTIR_objectFromFile',...
                        'AcquisitionData', arg);
                    self.read_only = true;
                end
            elseif isa(arg, 'sirf.STIR.AcquisitionData')
                self.handle_ = calllib...
                    ('mstir', 'mSTIR_acquisitionDataFromTemplate',...
                    arg.handle_);
            else
                error('AcquisitionData:wrong_ctor_source', ...
                'wrong source in AcquisitionData constructor')
            end
            sirf.Utilities.check_status(self.name, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function read_from_file(self, filename)
%***SIRF*** Reads acquisition data from a file.
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
            end
            self.handle_ = calllib('mstir', 'mSTIR_objectFromFile', ...
                'AcquisitionData', filename);
            sirf.Utilities.check_status(self.name, self.handle_);
            self.read_only = true;
        end
        function image = create_uniform_image(self, value, nx, ny)
%***SIRF*** create_uniform_image(value) creates compatible ImageData object.
%           The created object contains PET image of dimensions and voxel
%           sizes compatible with the scanner or its model used to produce
%           the acquisition data in self.
%           The specified value, if present, is assigned at all image voxels.
%           value: a float.
            image = sirf.STIR.ImageData();
            if nargin < 3
                image.handle_ = calllib...
                    ('mstir', 'mSTIR_imageFromAcquisitionData', self.handle_);
            else
                image.handle_ = calllib...
                    ('mstir', 'mSTIR_imageFromAcquisitionDataAndNxNy',...
                    self.handle_, nx, ny);
            end
            sirf.Utilities.check_status...
                ([self.name ':create_uniform_image'], image.handle_);
            if nargin > 1
                image.fill(value)
            end
        end
        function dim = dimensions(self)
%***SIRF*** Returns array of the acquisition data dimensions.
%           Dimensions are:
%           - number of tangential positions
%           - number of views
%           - number of (non-TOF) sinograms
%           - number of TOF bins
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib('mstir', 'mSTIR_getAcquisitionDataDimensions', ...
                self.handle_, ptr_i);
            dim = ptr_i.Value;
        end
        function data = as_array(self)
%***SIRF*** Returns 4D array of the acquisition data values.
%           Dimensions are:
%           - number of tangential positions
%           - number of views
%           - number of sinograms
%           - number of TOF bins
            dim = self.dimensions();
            n = dim(1)*dim(2)*dim(3)*dim(4);
            ptr_v = libpointer('singlePtr', zeros(n, 1));
            calllib('mstir', 'mSTIR_getAcquisitionData', self.handle_, ptr_v);
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3), dim(4));
        end
        function fill(self, value)
%***SIRF*** fill(value) fills the object with values;
%           value: float or array of floats (of the same dimensions and 
%                data order as the one returned by as_array() method) or an 
%                AcquisitionData object.
            if isempty(self.handle_)
                error([self.name ':fill'], ...
                    'AcquisitionData object not initialized')
            elseif self.read_only
                error([self.name ':fill'], ...
                    'Cannot fill read-only object, consider filling a clone')
            elseif isa(value, 'single')
                if numel(value) > 1
                    ptr_v = libpointer('singlePtr', value);
                    h = calllib('mstir', 'mSTIR_setAcquisitionData', ...
                        self.handle_, ptr_v);
                else
                    h = calllib('mstir', 'mSTIR_fillAcquisitionData', ...
                        self.handle_, value);
                end
                sirf.Utilities.check_status...
                    ([self.name ':fill'], h);
                sirf.Utilities.delete(h)
            elseif isa(value, 'double')
                if numel(value) > 1
                    ptr_v = libpointer('singlePtr', single(value));
                    h = calllib('mstir', 'mSTIR_setAcquisitionData', ...
                        self.handle_, ptr_v);
                else
                    h = calllib('mstir', 'mSTIR_fillAcquisitionData', ...
                        self.handle_, single(value));
                end
                sirf.Utilities.check_status...
                    ([self.name ':fill'], h);
                sirf.Utilities.delete(h)
            elseif isa(value, 'sirf.STIR.AcquisitionData')
                h = calllib('mstir', ...
                    'mSTIR_fillAcquisitionDataFromAcquisitionData', ...
                    self.handle_, value.handle_);
                sirf.Utilities.check_status([self.name ':fill'], h);
                sirf.Utilities.delete(h)
            else
                error([self.name ':fill'], 'wrong fill value')
            end
        end
        function ad = get_uniform_copy(self, value)
%***SIRF*** get_uniform_copy(value) returns a true copy of this object 
%           (not Matlab handle) filled with the specified value;
%           value: Matlab float
            if nargin < 2
                value = 0;
            end
            ad = sirf.STIR.AcquisitionData(self);
            ad.fill(value)
        end
    end
end
