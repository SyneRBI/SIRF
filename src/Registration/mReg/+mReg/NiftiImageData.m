classdef NiftiImageData < handle
% Class for image data.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018-2019 University College London
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
            name = 'NiftiImageData';
        end
    end
    methods
        function self = NiftiImageData(src)
            narginchk(0,1)
            self.name = 'NiftiImageData';
            if nargin < 1
                self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('mreg', 'mReg_objectFromFile', self.name, src);
            else
                error('NiftiImageData accepts no args or filename.')
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
            if isa(other, 'mReg.NiftiImageData')
            	calllib('mreg', 'mReg_NiftiImageData_maths_im', z.handle_, self.handle_, other.handle_, 0);
            elseif isnumeric(other)
            	calllib('mreg', 'mReg_NiftiImageData_maths_num', z.handle_, self.handle_, other, 0);
            else
            	error('NiftiImageData:plus should be image or number.')
            end
        	mUtilities.check_status('NiftiImageData:plus', z.handle_);
        end
        function z = minus(self, other)
            % Overloads the subtraction operator
            z = self.deep_copy();
            if isa(other, 'mReg.NiftiImageData')
            	calllib('mreg', 'mReg_NiftiImageData_maths_im', z.handle_, self.handle_, other.handle_, 1);
            elseif isnumeric(other)
            	calllib('mreg', 'mReg_NiftiImageData_maths_num', z.handle_, self.handle_, other, 1);
            else
            	error('NiftiImageData:minus should be image or number.')
            end
            mUtilities.check_status('NiftiImageData:minus', z.handle_);
        end
        function z = mtimes(self,other)
            % Overloads the multiplication operator
            if isnumeric(other)
            	z = self.deep_copy();
            	calllib('mreg', 'mReg_NiftiImageData_maths_num', z.handle_, self.handle_, other, 2);
            else
            	error('NiftiImageData:mtimes should be number.')
            end
            mUtilities.check_status('NiftiImageData:mtimes', z.handle_);
        end
        function value = eq(self, other)
        	% Overload equality operator
        	assert(isa(other, 'mReg.NiftiImageData'));
    	    hv = calllib('mreg', 'mReg_NiftiImageData_equal', self.handle_, other.handle_);
			mUtilities.check_status('parameter', hv);
    		value = logical(calllib('miutilities', 'mIntDataFromHandle', hv));
        end
        function value = ne(self, other)
        	% Overload inequality operator
        	value = ~(self==other);
        end
        function write(self, filename, datatype)
            %Save to file. See nifti1.h for datatypes (e.g., float (NIFTI_TYPE_FLOAT32) = 16)
            % Image's original datatpye is used by default.
            if nargin < 3
                datatype = -1;
            end
            h = calllib('mreg', 'mReg_NiftiImageData_write', self.handle_, filename, datatype);
            mUtilities.check_status([self.name ':write'], h);
            mUtilities.delete(h)
        end
        function value = get_max(self)
            %Get max.
            value = mReg.parameter(self.handle_, 'NiftiImageData', 'max', 'f');
        end
        function value = get_min(self)
            %Get min.
            value = mReg.parameter(self.handle_, 'NiftiImageData', 'min', 'f');
        end
        function value = get_sum(self)
            %Get sum.
            value = mReg.parameter(self.handle_, 'NiftiImageData', 'sum', 'f');
        end
        function value = get_dimensions(self)
            %Get dimensions.
            ptr_i = libpointer('int32Ptr', zeros(1, 8));
            calllib('mreg', 'mReg_NiftiImageData_get_dimensions', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function fill(self, val)
            %Fill image with single value.
            if numel(val) == 1
                h = calllib('mreg', 'mReg_NiftiImageData_fill', self.handle_, single(val));
            else
                if isa(val, 'single')
                    ptr_v = libpointer('singlePtr', val);
                else
                    ptr_v = libpointer('singlePtr', single(val));
                end
                assert(all(size(val) == size(self.as_array)),[self.name ':fill. Dimensions do not match.'])
                h = calllib('mreg', 'mReg_NiftiImageData_fill_arr', self.handle_, ptr_v);
            end
            mUtilities.check_status([self.name ':fill'], h);
            mUtilities.delete(h)            
        end
        function value = get_norm(self,other)
            %Get norm.
        	mUtilities.assert_validities(self,other)
    	    hv = calllib('mreg', 'mReg_NiftiImageData_norm', self.handle_, other.handle_);
			mUtilities.check_status('parameter', hv)
    		value = calllib('miutilities', 'mFloatDataFromHandle', hv);
        end
        function output = deep_copy(self)
            %Deep copy image.
            if strcmp(self.name, 'NiftiImageData')
                output = mReg.NiftiImageData();
            elseif strcmp(self.name, 'NiftiImageData3D')
                output = mReg.NiftiImageData3D();
            elseif strcmp(self.name, 'NiftiImageData3DTensor')
                output = mReg.NiftiImageData3DTensor();
            elseif strcmp(self.name, 'NiftiImageData3DDeformation')
                output = mReg.NiftiImageData3DDeformation();
            elseif strcmp(self.name, 'NiftiImageData3DDisplacement')
                output = mReg.NiftiImageData3DDisplacement();
            end
            calllib('mreg', 'mReg_NiftiImageData_deep_copy', output.handle_, self.handle_);
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function array = as_array(self)
            %Get data as numpy array.
            dim = self.get_dimensions();
            dim = dim(2:dim(1)+1);
            ptr_v = libpointer('singlePtr', zeros(dim));
            calllib('mreg', 'mReg_NiftiImageData_get_data', self.handle_, ptr_v);
            array = reshape(ptr_v.Value,dim);
        end
        function datatype = get_original_datatype(self)
            %Get original image datatype (internally everything is converted to float).
            h = calllib('mreg', 'mReg_NiftiImageData_get_original_datatype', self.handle_);
            mUtilities.check_status('NiftiImageData', h);
            datatype = calllib('miutilities', 'mCharDataFromHandle', h);
            mUtilities.delete(h)
        end
        function crop(self, min_, max_)
            assert(all(size(min_) == [1 7]), 'Min bounds should be a 1x7 array')
            assert(all(size(max_) == [1 7]), 'Max bounds should be a 1x7 array')
            min_ptr = libpointer('int32Ptr', single(min_));
            max_ptr = libpointer('int32Ptr', single(max_));
            h = calllib('mreg', 'mReg_NiftiImageData_crop', self.handle_, min_ptr, max_ptr);
            mUtilities.check_status('parameter', h)
        end
        function print_header(self)
            %Print metadata of nifti image.
            h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 1, self.handle_, [], [], [], []);
            mUtilities.check_status('parameter', h)
        end
    end
    methods(Static)
        function print_headers(to_print)
            %Print metadata of one or multiple (up to 5) nifti images.
            assert(ismatrix(to_print) && isa(to_print, 'mReg.NiftiImageData'), 'NiftiImageData.print_headers: give list of NiftiImageData.')
            num_ims = size(to_print,2);
            if num_ims == 1
                h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 1, to_print(1).handle_, [], [], [], []);
            elseif num_ims == 2
                h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 2, to_print(1).handle_, to_print(2).handle_, [], [], []);
            elseif num_ims == 3
                h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 3, to_print(1).handle_, to_print(2).handle_, to_print(3).handle_, [], []);
            elseif num_ims == 4
                h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 4, to_print(1).handle_, to_print(2).handle_, to_print(3).handle_, to_print(4).handle_, []);
            elseif num_ims == 5
                h = calllib('mreg', 'mReg_NiftiImageData_print_headers', 5, to_print(1).handle_, to_print(2).handle_, to_print(3).handle_, to_print(4).handle_, to_print(5).handle_);
            else
                error('print_headers only implemented for up to 5 images.')
            end
            mUtilities.check_status('parameter', h)
        end
    end
    % If you put this in, the workspace in matlab shows the size (eg., 64x64x64 NiftiImageData)
    % Without it, jusrt 1x1 NiftiImageData. However, with it, can't tell if we have an array of
    % NiftiImageData or if there is just one. TODO
    %methods (Hidden = true)
    %     function dim = size(self)
    %         % size
    %         dim = self.get_dimensions();
    %         dim = dim(2:dim(1)+1);
    %     end
    % end
end