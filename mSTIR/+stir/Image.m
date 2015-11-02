classdef Image < handle
    properties
        name
        handle
        rimsize
    end
    methods
        function self = Image(filename)
            self.name = 'Image';
            if nargin < 1
                error('Image:wrong_constructor', 'Wrong constructor for Image')
                %self.handle = [];
            else
                self.handle = calllib('mstir', 'mSTIR_imageFromFile', filename);
            end
            self.rimsize = -1;
        end
        function delete(self)
            %if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle, self.name)
            %end
        end
%         function read_from_file(self, filename)
%             if ~isempty(self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle, self.name)
%             end
%             self.handle = calllib('mstir', 'mSTIR_imageFromFile', filename);
%             stir.checkExecutionStatus('read_from_file', self.handle);
%         end
        function diff = diff_from(self, image)
            h = calllib('mstir', 'mSTIR_imagesDifference',...
                     self.handle, image.handle, self.rimsize);
            stir.checkExecutionStatus('diff_from', h);
            diff = calllib('mstir', 'mDoubleDataFromHandle', h);
            calllib('mstir', 'mDeleteDataHandle', h)
        end
        function data = density(self)
            [ptr, dim] = calllib...
                ('mstir', 'mSTIR_getImageDimensions', self.handle, zeros(3, 1));
            n = dim(1)*dim(2)*dim(3);
            [ptr, data] = calllib...
                ('mstir', 'mSTIR_getImageData', self.handle, zeros(n, 1));
            data = reshape(data, dim(3), dim(2), dim(1));
        end
    end
end