classdef ImageData < handle
% Class for PET image data objects.

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
        handle
        rimsize
    end
    methods
        function self = ImageData(filename)
            % Creates an ImageData object. 
            % The optional argument is the name of an Interfile containing
            % image data.
            % If no argument given, the object remains empty, and needs to
            % be defined by its initialise() method before it can be used.
            self.name = 'ImageData';
            if nargin < 1
                self.handle = [];
            else
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
                mUtil.checkExecutionStatus('ImageData', self.handle)
            end
            self.rimsize = -1;
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function initialise(self,...
                arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
%***SIRF*** Sets this image size in voxels, voxel sizes in mm and the origin.
%         All arguments except the first one are optional.
%         Present arguments are either all scalars or all 3-component arrays.
%         The first array argument or three scalar arguments set the image
%         sizes in voxels.
%         The second array argument or three scalar arguments set the voxel
%         sizes in mm (if absent, sizes default to (1,1,1)).
%         The third array argument or three scalar arguments set the origin
%         (if absent, defaults to (0,0,0)).
            vsize = [1 1 1];
            origin = [0 0 0];
            if max(size(arg1)) == 1
                dim = [arg1 arg2 arg3];
                if nargin > 4
                    vsize = [arg4 arg5 arg6];
                    if nargin > 7
                        origin = [arg7, arg8, arg9];
                    end
                end
            else
                dim = arg1;
                if nargin > 2
                    vsize = arg2;
                    if nargin > 3
                        origin = arg3;
                    end
                end
            end
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
            voxels = calllib('mstir', 'mSTIR_voxels3DF',...
                dim(1), dim(2), dim(3),...
                vsize(1), vsize(2), vsize(3),...
                origin(1), origin(2), origin(3));
            mUtil.checkExecutionStatus('ImageData:initialise', voxels)
            self.handle = calllib('mstir', 'mSTIR_imageFromVoxels', voxels);
            mUtil.checkExecutionStatus('ImageData:initialise', self.handle)
            calllib('mutilities', 'mDeleteDataHandle', voxels)
        end
        function fill(self, value)
%***SIRF*** Sets this image values at voxels.
%         The argument is either 3D array of values or a scalar to be
%         assigned at each voxel.
            if numel(value) == 1
                calllib('mstir', 'mSTIR_fillImage', self.handle, value)
            else
                ptr_v = libpointer('doublePtr', value);
                calllib('mstir', 'mSTIR_setImageData', self.handle, ptr_v)
            end
        end
        function image = clone(self)
%***SIRF*** Creates a copy of this image.
            image = mStir.ImageData();
            image.handle = calllib('mstir', 'mSTIR_imageFromImage',...
                self.handle);
            mUtil.checkExecutionStatus('ImageData:clone', self.handle)
        end
        function image = get_uniform_copy(self, value)
%***SIRF*** Creates a copy of this image filled with the specified value.
            if nargin < 2
                value = 1.0;
            end
            image = mStir.ImageData();
            image.handle = calllib('mstir', 'mSTIR_imageFromImage',...
                self.handle);
            mUtil.checkExecutionStatus('ImageData:get_uniform_copy', self.handle)
            image.fill(value)
        end
        function read_from_file(self, filename)
%***SIRF*** Reads the image data from a file.
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
            self.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
            mUtil.checkExecutionStatus('ImageData:read_from_file', self.handle);
        end
        function write(self, filename)
            h = calllib('mstir', 'mSTIR_writeImage', self.handle, filename);
            mUtil.checkExecutionStatus('ImageData:write', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function add_shape(self, shape, add)
%***SIRF*** Adds a uniform shape to the image. 
%         The image values at voxels inside the added shape are increased 
%         by the value of the last argument.
            if isempty(self.handle)
                error('ImageData:error', 'cannot add shapes to uninitialised image');
            end
            h = calllib...
                ('mstir', 'mSTIR_addShape', self.handle,...
                shape.handle, add);
            mUtil.checkExecutionStatus('ImageData:add_shape', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function diff = diff_from(self, image)
%***SIRF*** Returns the relative difference between self and the image
%         specified by the last argument, i.e. the maximal difference at
%         voxels of common containing box divided by the maximum value
%         of self.
            h = calllib('mstir', 'mSTIR_imagesDifference',...
                     self.handle, image.handle, self.rimsize);
            mUtil.checkExecutionStatus('ImageData:diff_from', h);
            diff = calllib('mutilities', 'mDoubleDataFromHandle', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function data = as_array(self)
%***SIRF*** Returns 3D array of this image values at voxels.

%             [ptr, dim] = calllib...
%                 ('mstir', 'mSTIR_getImageDimensions', self.handle, zeros(3, 1));
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib...
                ('mstir', 'mSTIR_getImageDimensions', self.handle, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
%             [ptr, data] = calllib...
%                 ('mstir', 'mSTIR_getImageData', self.handle, zeros(n, 1));
%             data = reshape(data, dim(3), dim(2), dim(1));
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mstir', 'mSTIR_getImageData', self.handle, ptr_v)
            data = reshape(ptr_v.Value, dim(3), dim(2), dim(1));
        end
        function show(self)
%***SIRF*** Interactively plots this image data as a set of 2D image slices.
            data = self.as_array();
            shape = size(data);
            nz = shape(3);
            if nz < 1
                return
            end
            data = data/max(data(:));
            fprintf('Please enter the number of the image slice to view\n')
            fprintf('(a value outside the range [1 : %d] will stop this loop)\n', nz)
            while true
                z = input('image: ');
                if z < 1 || z > nz
                    break
                end
                figure(z)
                imshow(data(:, :, z))
                title(['image ' num2str(z)])
            end
        end
    end
end