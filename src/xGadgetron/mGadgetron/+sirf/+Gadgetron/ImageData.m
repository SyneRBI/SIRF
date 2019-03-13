classdef ImageData < mSIRF.ImageData
% Class for MR image data objects.

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
        name_
    end
    methods (Static)
        function name = class_name()
            name = 'ImageData';
        end
        function obj = same_object()
            obj = sirf.Gadgetron.ImageData();
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
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function read_from_file(self, file)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
            end
            self.handle_ = calllib('mgadgetron', 'mGT_readImages', file);
            mUtilities.check_status(self.name_, self.handle_);
        end
        function img = image(self, num)
            img = sirf.Gadgetron.Image();
            img.handle_ = calllib('mgadgetron', 'mGT_imageWrapFromContainer', ...
                self.handle_, num - 1);
            mUtilities.check_status(self.name_, img.handle_);
        end
        function images = select(self, attribute, value)
%***STIR*** select(attribute, value) returns a subset of this image data 
%         with the specified value of the specified attribute;
%         attribute: attribute name (Matlab char string)
%         value    : attribute value (Matlab char string)
            if isempty(self.handle_)
                error('ImageData:empty_object', ...
                    'cannot handle empty object')
            end
            images = sirf.Gadgetron.ImageData();
            images.handle_ = calllib('mgadgetron', 'mGT_selectImages', ...
                self.handle_, attribute, value);
            mUtilities.check_status(self.name_, images.handle_);
        end
        function images = process(self, list)
%***SIRF*** process(list) returns images processed by a chain of gadgets.
%             The argument is a cell array of gadget definitions
%             [{gadget1_definition}, {gadget2_definition}, ...],
%             where gadget definitions are strings of the form
%             'label:name(property1=value1,property2=value2,...)',
%             where the only mandatory field is name, the Gadgetron 
%             name of the gadget. An optional expression in round 
%             brackets can be used to assign values to gadget properties,
%             and an optional label can be used to change the labelled
%             gadget properties after the chain has been defined.
            ip = sirf.Gadgetron.ImageDataProcessor(list);
            images = ip.process(self);
        end
        function ft = is_real(self)
%***SIRF*** Returns true if this image data is real and false otherwise.
            handle = calllib('mgadgetron', 'mGT_imageDataType', ...
                self.handle_, 0);
            mUtilities.check_status(self.name_, handle);
            v = calllib('miutilities', 'mIntDataFromHandle', handle);
            %v = calllib('mutilities', 'mIntDataFromHandle', handle);
            %calllib('mutilities', 'mDeleteDataHandle', handle)
            mUtilities.delete(handle)
            ft = (v ~= 7 && v ~= 8);
        end
        function data = as_array(self)
%***SIRF*** Returns 3D complex array representing this image data.
%         First two dimensions are x and y, the third is a product of all
%         other dimensions (z/slice/repetition etc.).
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            if self.number() > 0
                img = self.image(1);
                calllib('mgadgetron', 'mGT_getImageDim', img.handle_, ptr_i);
            else
                data = [];
                return
            end
            dim = ptr_i.Value;
            nz = dim(3)*dim(4)*self.number();
            n = dim(1)*dim(2)*nz;
            if self.is_real()
                ptr_v = libpointer('singlePtr', zeros(n, 1));
                h = calllib('mgadgetron', 'mGT_getImagesDataAsFloatArray', ...
                    self.handle_, ptr_v);
                data = reshape(ptr_v.Value, dim(1), dim(2), nz);
            else
                ptr_z = libpointer('singlePtr', zeros(2, n));
                h = calllib...
                    ('mgadgetron', 'mGT_getImagesDataAsCmplxArray', ...
                    self.handle_, ptr_z);
                data = reshape(ptr_z.Value(1:2:end) + 1i*ptr_z.Value(2:2:end), ...
                    dim(1), dim(2), nz);
            end
            mUtilities.check_status('ImageData', h);
            mUtilities.delete(h)
        end
        function fill(self, data)
%***SIRF*** Changes image data to that in 3D complex array argument.
            if isempty(self.handle_)
                error('ImageData:empty_object', 'cannot handle empty object')
            end
            if self.is_real()
                h = calllib('mgadgetron', 'mGT_setImagesDataAsRealArray', ...
                    self.handle_, data);
            else
                z = [real(data(:))'; imag(data(:))'];
                if isa(z, 'single')
                    ptr_z = libpointer('singlePtr', z);
                else
                    ptr_z = libpointer('singlePtr', single(z));
                end
                h = calllib('mgadgetron', 'mGT_setImagesDataAsCmplxArray', ...
                    self.handle_, ptr_z);
            end
            mUtilities.check_status('ImageData', h);
            mUtilities.delete(h)
            %calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function show(self)
%***SIRF*** Interactively plots this image data as a set of 2D image slices.
            nz = self.number();
            if nz < 1
                return
            end
            data = self.as_array();
            if ~self.is_real()
                data = abs(data);
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