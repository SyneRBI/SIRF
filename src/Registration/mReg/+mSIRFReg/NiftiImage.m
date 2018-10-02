classdef NiftiImage < handle
% Class for image data.

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
        handle_
    end
    methods(Static)
        function name = class_name()
            name = 'NiftiImage';
        end
    end
    methods
        function self = NiftiImage(src)
            narginchk(0,1)
            self.name = 'NiftiImage';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, src);
            else
                error('NiftiImage accepts no args or filename.')
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function z = plus(self, other)
            % Overloads the addition operator
            z = self.deep_copy();
            if isa(other, 'mSIRFReg.NiftiImage')
            	calllib('msirfreg', 'mSIRFReg_NiftiImage_maths_im', z.handle_, self.handle_, other.handle_, 0);
            elseif isnumeric(other)
            	calllib('msirfreg', 'mSIRFReg_NiftiImage_maths_num', z.handle_, self.handle_, other, 0);
            else
            	error('NiftiImage:plus should be image or number.')
            end
        	mUtilities.check_status('NiftiImage:plus', z.handle_);
        end
        function z = minus(self, other)
            % Overloads the subtraction operator
            z = self.deep_copy();
            if isa(other, 'mSIRFReg.NiftiImage')
            	calllib('msirfreg', 'mSIRFReg_NiftiImage_maths_im', z.handle_, self.handle_, other.handle_, 1);
            elseif isnumeric(other)
            	calllib('msirfreg', 'mSIRFReg_NiftiImage_maths_num', z.handle_, self.handle_, other, 1);
            else
            	error('NiftiImage:minus should be image or number.')
            end
            mUtilities.check_status('NiftiImage:minus', z.handle_);
        end
        function z = mtimes(self,other)
            % Overloads the multiplication operator
            if isnumeric(other)
            	z = self.deep_copy();
            	calllib('msirfreg', 'mSIRFReg_NiftiImage_maths_num', z.handle_, self.handle_, other, 2);
            else
            	error('NiftiImage:mtimes should be number.')
            end
            mUtilities.check_status('NiftiImage:mtimes', z.handle_);
        end
        function value = eq(self, other)
        	% Overload equality operator
        	assert(isa(other, 'mSIRFReg.NiftiImage'));
    	    hv = calllib('msirfreg', 'mSIRFReg_NiftiImage_equal', self.handle_, other.handle_);
			mUtilities.check_status('parameter', hv);
    		value = logical(calllib('miutilities', 'mIntDataFromHandle', hv));
        end
        function value = ne(self, other)
        	% Overload inequality operator
        	value = ~(self==other);
        end
        function save_to_file(self, filename, datatype)
            %Save to file.
            if nargin < 3
                datatype ='';
            end
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage_save_to_file', self.handle_, filename, datatype);
            mUtilities.check_status([self.name ':save_to_file'], h);
            mUtilities.delete(h)
        end
        function value = get_max(self)
            %Get max.
            value = mSIRFReg.parameter(self.handle_, 'NiftiImage', 'max', 'f');
        end
        function value = get_min(self)
            %Get min.
            value = mSIRFReg.parameter(self.handle_, 'NiftiImage', 'min', 'f');
        end
        function value = get_sum(self)
            %Get sum.
            value = mSIRFReg.parameter(self.handle_, 'NiftiImage', 'sum', 'f');
        end
        function value = get_dimensions(self)
            %Get dimensions.
            ptr_i = libpointer('int32Ptr', zeros(1, 8));
            calllib('msirfreg', 'mSIRFReg_NiftiImage_get_dimensions', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function fill(self, val)
            %Fill image with single value.
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage_fill', self.handle_, val);
            mUtilities.check_status([self.name ':fill'], h);
            mUtilities.delete(h)            
        end
        function value = get_norm(self,other)
            %Get norm.
        	mUtilities.assert_validities(self,other)
    	    hv = calllib('msirfreg', 'mSIRFReg_NiftiImage_norm', self.handle_, other.handle_);
			mUtilities.check_status('parameter', hv)
    		value = calllib('miutilities', 'mFloatDataFromHandle', hv);
        end
        function output = deep_copy(self)
            %Deep copy image.
            if strcmp(self.name, 'NiftiImage')
                output = mSIRFReg.NiftiImage();
            elseif strcmp(self.name, 'NiftiImage3D')
                output = mSIRFReg.NiftiImage3D();
            elseif strcmp(self.name, 'NiftiImage3DTensor')
                output = mSIRFReg.NiftiImage3DTensor();
            elseif strcmp(self.name, 'NiftiImage3DDeformation')
                output = mSIRFReg.NiftiImage3DDeformation();
            elseif strcmp(self.name, 'NiftiImage3DDisplacement')
                output = mSIRFReg.NiftiImage3DDisplacement();
            end
            calllib('msirfreg', 'mSIRFReg_NiftiImage_deep_copy', output.handle_, self.handle_);
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function array = as_array(self)
            %Get data as numpy array.
            dim = self.get_dimensions();
            dim = dim(2:dim(1)+1);
            ptr_v = libpointer('singlePtr', zeros(dim));
            calllib('msirfreg', 'mSIRFReg_NiftiImage_get_data', self.handle_, ptr_v);
            array = reshape(ptr_v.Value,dim);
        end
        function datatype = get_original_datatype(self)
            %Get original image datatype (internally everything is converted to float).
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage_get_original_datatype', self.handle_);
            mUtilities.check_status('NiftiImage', h);
            datatype = calllib('miutilities', 'mCharDataFromHandle', h);
            mUtilities.delete(h)
        end
        function crop(self, min_, max_)
            assert(all(size(min_) == [1 7]), 'Min bounds should be a 1x7 array')
            assert(all(size(max_) == [1 7]), 'Max bounds should be a 1x7 array')
            min_ptr = libpointer('int32Ptr', single(min_));
            max_ptr = libpointer('int32Ptr', single(max_));
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage_crop', self.handle_, min_ptr, max_ptr);
            mUtilities.check_status('parameter', h)
        end
        function dump_header(self)
            %Dump metadata of nifti image.
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 1, self.handle_, [], [], [], []);
            mUtilities.check_status('parameter', h)
        end
    end
    methods(Static)
        function dump_headers(to_dump)
            %Dump metadata of one or multiple (up to 5) nifti images.
            assert(ismatrix(to_dump) && isa(to_dump, 'mSIRFReg.NiftiImage'), 'NiftiImage.dump_headers: give list of NiftiImage.')
            num_ims = size(to_dump,2);
            if num_ims == 1
                h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 1, to_dump(1).handle_, [], [], [], []);
            elseif num_ims == 2
                h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 2, to_dump(1).handle_, to_dump(2).handle_, [], [], []);
            elseif num_ims == 3
                h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 3, to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_, [], []);
            elseif num_ims == 4
                h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 4, to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_, to_dump(4).handle_, []);
            elseif num_ims == 5
                h = calllib('msirfreg', 'mSIRFReg_NiftiImage_dump_headers', 5, to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_, to_dump(4).handle_, to_dump(5).handle_);
            else
                error('dump_nifti_info only implemented for up to 5 images.')
            end
            mUtilities.check_status('parameter', h)
        end
    end
    % If you put this in, the workspace in matlab shows the size (eg., 64x64x64 NiftiImage)
    % Without it, jusrt 1x1 NiftiImage. However, with it, can't tell if we have an array of
    % NiftiImage or if there is just one. TODO
    %methods (Hidden = true)
    %     function dim = size(self)
    %         % size
    %         dim = self.get_dimensions();
    %         dim = dim(2:dim(1)+1);
    %     end
    % end
end