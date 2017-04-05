classdef ImageData < mGadgetron.DataContainer
% Class for MR image data objects.

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
        name_
    end
    methods (Static)
        function obj = same_object()
            obj = mGadgetron.ImageData();
        end
    end
    methods
        function self = ImageData()
            % Creates empty ImageData object.
            self.name_ = 'ImageData';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function write(self, file, group)
%***SIRF*** write(file, group) writes this image to a file in HDF5 format;
%         file : file name (Matlab char string)
%         group: group name (Matlab char string)
            if isempty(self.handle_)
                error('ImageData:empty_object', ...
                    'cannot handle empty object')
            end
            handle = calllib('mgadgetron', 'mGT_writeImages', ...
                self.handle_, file, group);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function images = select(self, attr, value)
            % Returns images with given value of given attribute.
            if isempty(self.handle_)
                error('ImageData:empty_object', ...
                    'cannot handle empty object')
            end
            images = mGadgetron.ImageData();
            images.handle_ = calllib('mgadgetron', 'mGT_selectImages', ...
                self.handle_, attr, value);
            mUtil.checkExecutionStatus(self.name_, images.handle_);
        end
        function images = process(self, list)
            % Returns images processed by a chain of gadgets.
            % The argument is a cell array of gadget definitions
            % [{gadget1_definition}, {gadget2_definition}, ...],
            % where gadget definitions are strings of the form
            % 'label:name(property1=value1,property2=value2,...)',
            % where the only mandatory field is name, the Gadgetron 
            % name of the gadget. An optional expression in round 
            % brackets can be used to assign values to gadget properties,
            % and an optional label can be used to change the labelled
            % gadget properties after the chain has been defined.
            ip = mGadgetron.ImageDataProcessor(list);
            images = ip.process(self);
        end
        function ft = is_real(self)
            handle = calllib('mgadgetron', 'mGT_imageDataType', ...
                self.handle_, 0);
            mUtil.checkExecutionStatus(self.name_, handle);
            v = calllib('mutilities', 'mIntDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            ft = (v ~= 7 && v ~= 8);
        end
        function data = as_array(self)
            % Returns 3D complex array representing 3D image.
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            if self.number() > 0
                calllib('mgadgetron', 'mGT_getImageDimensions', ...
                    self.handle_, 0, ptr_i);
            else
                data = [];
                return
            end
            dim = ptr_i.Value;
            nz = dim(3)*dim(4)*self.number();
            n = dim(1)*dim(2)*nz;
            if self.is_real()
                ptr_v = libpointer('doublePtr', zeros(n, 1));
                calllib('mgadgetron', 'mGT_getImagesDataAsDoubleArray', ...
                    self.handle_, ptr_v)
                data = reshape(ptr_v.Value, dim(1), dim(2), nz);
            else
                ptr_re = libpointer('doublePtr', zeros(n, 1));
                ptr_im = libpointer('doublePtr', zeros(n, 1));
                calllib...
                    ('mgadgetron', 'mGT_getImagesDataAsComplexArray', ...
                    self.handle_, ptr_re, ptr_im)
                re = reshape(ptr_re.Value, dim(1), dim(2), nz);
                im = reshape(ptr_im.Value, dim(1), dim(2), nz);
                data = re + 1j*im;
            end
        end
        function show(self)
            % Plots 2D images (slices).
            ni = self.number();
            if ni < 1
                return
            end
            data = self.as_array();
            if ~self.is_real()
                data = abs(data);
            end
            data = data/max(max(max(data)));
            fprintf('Please enter the number of the image to view\n')
            fprintf('(a value outside the range [1 : %d] will stop this loop)\n', ni)
            while true
                i = input('image: ');
                if i < 1 || i > ni
                    break
                end
                figure(i)
                imshow(data(:, :, i))
                title(['image ' num2str(i)])
            end
        end
    end
end