classdef ImageData < sirf.SIRF.ImageData
% Class for PET image data objects.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2020 University College London.
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
        name
        rimsize
    end
    methods(Static)
        function name = class_name()
            name = 'ImageData';
        end
        function obj = same_object()
            obj = sirf.STIR.ImageData();
        end
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
                self.handle_ = [];
            else
                self.handle_ = calllib...
                    ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
                sirf.Utilities.check_status('ImageData', self.handle_)
            end
            self.rimsize = -1;
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
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
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
            end
            voxels = calllib('mstir', 'mSTIR_voxels3DF',...
                dim(1), dim(2), dim(3),...
                vsize(1), vsize(2), vsize(3),...
                origin(1), origin(2), origin(3));
            sirf.Utilities.check_status('ImageData:initialise', voxels)
            self.handle_ = calllib('mstir', 'mSTIR_imageFromVoxels', voxels);
            sirf.Utilities.check_status('ImageData:initialise', self.handle_)
            sirf.Utilities.delete(voxels)
        end
        function fill(self, value)
%***SIRF*** Sets this image values at voxels.
%         The argument is either 3D array of values or a scalar to be
%         assigned at each voxel.
            if numel(value) == 1
                h = calllib('mstir', 'mSTIR_fillImage', ...
                    self.handle_, single(value));
            else
                if isa(value, 'single')
                    ptr_v = libpointer('singlePtr', value);
                else
                    ptr_v = libpointer('singlePtr', single(value));
                end
                h = calllib('mstir', 'mSTIR_setImageData', self.handle_, ptr_v);
            end
            sirf.Utilities.check_status('ImageData:fill', h)
            sirf.Utilities.delete(h)
        end
        function image = get_uniform_copy(self, value)
%***SIRF*** Creates a copy of this image filled with the specified value.
            if nargin < 2
                value = 1.0;
            end
            image = sirf.STIR.ImageData();
            image.handle_ = calllib('mstir', 'mSTIR_imageFromImage',...
                self.handle_);
            sirf.Utilities.check_status('ImageData:get_uniform_copy', self.handle_)
            image.fill(value)
        end
        function read_from_file(self, filename)
%***SIRF*** Reads the image data from a file.
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                sirf.Utilities.delete(self.handle_)
            end
            self.handle_ = calllib...
                ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
            sirf.Utilities.check_status('ImageData:read_from_file', self.handle_);
        end
        function add_shape(self, shape, add, num_samples_in_each_direction)
%***SIRF*** Adds a uniform shape to the image. 
%         The image values at voxels inside the added shape are increased 
%         by the value of the last argument.
%
%         If a shape partially fills a voxel, it is possible to choose the
%         number of samples that will be used in each direction to determine the
%         fraction of the voxel that is filled by the shape. For a 3D image,
%         using num_samples_in_each_direction=2 would result in 2^3=8 samples.
            sirf.Utilities.assert_validity(shape, 'Shape')
            if isempty(self.handle_)
                error('ImageData:error', 'cannot add shapes to uninitialised image');
            end
            if nargin < 4
                num_samples_in_each_direction = 1
            end
            h = calllib...
                ('mstir', 'mSTIR_addShape', self.handle_,...
                shape.handle_, add, num_samples_in_each_direction);
            sirf.Utilities.check_status('ImageData:add_shape', h);
            sirf.Utilities.delete(h)
        end
        function dim = size(self)
%***SIRF*** Returns the dimensions of 3D array of this image values at voxels.
            if isempty(self.handle_)
                dim = [];
                return
            end
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            h = calllib...
                ('mstir', 'mSTIR_getImageDimensions', self.handle_, ptr_i);
            sirf.Utilities.check_status('ImageData:as_array', h);
            sirf.Utilities.delete(h)
            idim = ptr_i.Value;
            dim = [idim(3) idim(2) idim(1)];
%             dim = fliplr(idim); % does not work
        end
        function data = as_array(self)
%***SIRF*** Returns 3D array of this image values at voxels.

%             [ptr, dim] = calllib...
%                 ('mstir', 'mSTIR_getImageDimensions', self.handle_, zeros(3, 1));
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            h = calllib...
                ('mstir', 'mSTIR_getImageDimensions', self.handle_, ptr_i);
            sirf.Utilities.check_status('ImageData:as_array', h);
            sirf.Utilities.delete(h)
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
%             [ptr, data] = calllib...
%                 ('mstir', 'mSTIR_getImageData', self.handle_, zeros(n, 1));
%             data = reshape(data, dim(3), dim(2), dim(1));
            ptr_v = libpointer('singlePtr', zeros(n, 1));
            h = calllib...
                ('mstir', 'mSTIR_getImageData', self.handle_, ptr_v);
            sirf.Utilities.check_status('ImageData:as_array', h);
            sirf.Utilities.delete(h)
            data = reshape(ptr_v.Value, dim(3), dim(2), dim(1));
        end
        function write(self,filename,par)
            if nargin < 3
                write@sirf.SIRF.ImageData(self, filename)
                return
            end
        %Write with parameter file
            h = calllib...
                ('mstir', 'mSTIR_writeImage_par', self.handle_, filename, par);
            sirf.Utilities.check_status('ImageData:write_w_param_file', h);
            sirf.Utilities.delete(h);
        end
        function show(self, z)
%***SIRF*** Interactively plots this image data as a set of 2D image slices.
            data = self.as_array();
            shape = size(data);
            nz = shape(3);
            if nz < 1
                return
            end
            if nargin > 1
                the_title = sprintf('Slice %d', z);
                sirf.Utilities.show_2D_array(data(:,:,z), the_title, 'x', 'y');
                return
            end                
            data = data/max(data(:));
            fprintf('Please enter z-slice numbers (ex: 1, 3-5) %s\n', ...
                'or 0 to stop the loop')
            while true
                s = input('z-slices to display: ', 's');
                err = sirf.Utilities.show_3D_array...
                    (data, 'Selected slices', 'x', 'y', 'slice', s); %elect);
                if err ~= 0
                    fprintf('out-of-range slice numbers selected, %s\n', ...
                        'quitting the loop')
                    break
                end
            end
        end
    end
end
