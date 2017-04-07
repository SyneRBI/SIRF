classdef AcquisitionData < handle
% Class for PET acquisition data objects.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
        handle
    end
    methods
        function self = AcquisitionData(arg)
%         AcquisitionData(arg) creates new AcquisitionData object 
%         from a file or another AcquisitionData object;
%         arg:  file name or AcquisitionData object.
            self.handle = [];
            self.name = 'AcquisitionData';
            if nargin < 1
                return
            elseif ischar(arg)
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile',...
                    'AcquisitionData', arg);
            elseif isa(arg, 'mStir.AcquisitionData')
                self.handle = calllib...
                    ('mstir', 'mSTIR_acquisitionsDataFromTemplate',...
                    arg.handle);
            else
                error('AcquisitionData:wrong_ctor_source', ...
                'wrong source in AcquisitionData constructor')
            end
            mUtil.checkExecutionStatus(self.name, self.handle);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function read_from_file(self, filename)
%***SIRF*** Reads acquisition data from a file.
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
            self.handle = calllib('mstir', 'mSTIR_objectFromFile', ...
                'AcquisitionData', filename);
            mUtil.checkExecutionStatus(self.name, self.handle);
        end
        function image = create_uniform_image(self, value)
%***SIRF*** create_uniform_image(value) creates compatible ImageData object.
%         The created object contains PET image of dimensions and voxel sizes 
%         compatible with the scanner or its model used to produce the
%         acquisition data in self.
%         The specified value, if present, is assigned at all image voxels.
%         value: a double.
            image = mStir.ImageData();
            image.handle = calllib...
                ('mstir', 'mSTIR_imageFromAcquisitionData', self.handle);
            mUtil.checkExecutionStatus...
                ([self.name ':create_uniform_image'], image.handle);
            if nargin > 1
                image.fill(value)
            end
        end
        function data = as_array(self)
%***SIRF*** Returns 3D array of the acquisition data values.
%         Dimensions are:
%         - number of tangential positions
%         - number of views
%         - number of sinograms
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsDimensions', ...
                self.handle, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsData', self.handle, ptr_v);
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3));
        end
        function fill(self, value)
%***SIRF*** fill(value) fills the object with values;
%         value: double or array of doubles (of the same dimensions and data
%                order as the one returned by as_array() method) or an 
%                AcquisitionData object.
            if isempty(self.handle)
                error([self.name ':fill'], ...
                    'AcquisitionData object not initialized')
            elseif isa(value, 'double')
                if numel(value) > 1
                    ptr_v = libpointer('doublePtr', value);
                    calllib('mstir', 'mSTIR_setAcquisitionsData', ...
                        self.handle, ptr_v);
                else
                    calllib('mstir', 'mSTIR_fillAcquisitionsData', ...
                        self.handle, value);
                end
            elseif isa(value, 'mStir.AcquisitionData')
                calllib('mstir', ...
                    'mSTIR_fillAcquisitionsDataFromAcquisitionsData', ...
                    self.handle, value.handle);
            else
                error([self.name ':fill'], 'wrong fill value')
            end
        end
        function ad = clone(self)
%***SIRF*** Returns a true copy of this object (not Matlab handle).
            ad = mStir.AcquisitionData(self);
            ad.fill(self)
        end
        function ad = get_uniform_copy(self, value)
%***SIRF*** get_uniform_copy(value) returns a true copy of this object 
%         (not Matlab handle) filled with the specified value;
%         value: Matlab double
            if nargin < 2
                value = 0;
            end
            ad = mStir.AcquisitionData(self);
            ad.fill(value)
        end
    end
end