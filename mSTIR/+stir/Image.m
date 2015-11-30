classdef Image < handle
    properties
        name
        handle
        voxels
        rimsize
    end
    methods
        function self = Image(filename)
            self.name = 'Image';
            if nargin < 1
                self.handle = [];
            else
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
%                self.handle = calllib('mstir', 'mSTIR_imageFromFile', filename);
                stir.checkExecutionStatus('Image', self.handle)
            end
            self.voxels = [];
            self.rimsize = -1;
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
            if ~isempty(self.voxels)
                calllib('mstir', 'mSTIR_deleteObject', self.voxels)
            end
        end
        function initialise(self,...
                arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
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
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
            if ~isempty(self.voxels)
                calllib('mstir', 'mSTIR_deleteObject', self.voxels)
            end
            self.voxels = calllib('mstir', 'mSTIR_voxels3DF',...
                dim(1), dim(2), dim(3),...
                vsize(1), vsize(2), vsize(3),...
                origin(1), origin(2), origin(3));
            stir.checkExecutionStatus('Image.initialise', self.voxels)
            self.handle = calllib('mstir', 'mSTIR_imageFromVoxels',...
                self.voxels);
            stir.checkExecutionStatus('Image.initialise', self.handle)
        end
        function fill(self, value)
            calllib('mstir', 'mSTIR_fillImage', self.handle, value)
        end
        function image = clone(self)
            image = stir.Image();
            image.handle = calllib('mstir', 'mSTIR_imageFromImage',...
                self.handle);
            stir.checkExecutionStatus('Image.clone', self.handle)
        end
        function image = get_empty_copy(self, value)
            if nargin < 2
                value = 1.0;
            end
            image = stir.Image();
            image.handle = calllib('mstir', 'mSTIR_imageFromImage',...
                self.handle);
            stir.checkExecutionStatus('Image.get_empty_copy', self.handle)
            image.fill(value)
        end
        function read_from_file(self, filename)
            if ~isempty(self.handle)
               calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
            if ~isempty(self.voxels)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile', 'Image', filename);
%            self.handle = calllib('mstir', 'mSTIR_imageFromFile', filename);
            stir.checkExecutionStatus('Image.read_from_file', self.handle);
        end
        function add_shape(self, shape, scale)
            if isempty(self.handle)
                error('Image:error', 'cannot add shapes to uninitialised image');
            end
            h = calllib...
                ('mstir', 'mSTIR_addShape', self.handle, self.voxels,...
                shape.handle, scale);
            stir.checkExecutionStatus('Image.add_shape', h);
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function diff = diff_from(self, image)
            h = calllib('mstir', 'mSTIR_imagesDifference',...
                     self.handle, image.handle, self.rimsize);
            stir.checkExecutionStatus('Image.diff_from', h);
            diff = calllib('mstir', 'mDoubleDataFromHandle', h);
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function data = density(self)
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
    end
end