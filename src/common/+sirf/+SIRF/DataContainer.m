classdef DataContainer < handle
% INTERNAL USE ONLY.
% Class for an abstract data container.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
% Copyright 2020 University College London
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
        handle_
    end
    methods (Abstract, Static)
        same_object()
    end
    methods
        function self = DataContainer()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function num = number(self)
%***SIRF*** Returns the number of items in this container.
            handle = calllib('msirf', 'mSIRF_dataItems', self.handle_);
            sirf.Utilities.check_status('DataContainer', handle);
            num = calllib('miutilities', 'mIntDataFromHandle', handle);
            sirf.Utilities.delete(handle)
        end
        function empty = is_empty(self)
%***SIRF*** Returns true if this container does not actually store data
%           (but may store metadata - geometry etc.), e.g. is a template.
            empty = self.number() < 1;
        end
        function copy = clone(self)
            if isempty(self.handle_)
                error('DataContainer:clone:empty_object', ...
                    'cannot handle empty object')
            end
            copy = self.same_object();
            copy.handle_ = calllib('msirf', 'mSIRF_clone', self.handle_);
            sirf.Utilities.check_status('DataContainer', copy.handle_);
        end
        function r = norm(self)
%***SIRF*** Returns the 2-norm of this data container viewed as a vector.
            handle = calllib('msirf', 'mSIRF_norm', self.handle_);
            sirf.Utilities.check_status('DataContainer', handle);
            r = calllib('miutilities', 'mFloatDataFromHandle', handle);
            sirf.Utilities.delete(handle)
        end
        function z = dot(self, other)
%***SIRF*** Returns the dot product of this data container with another one 
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            handle = calllib('msirf', 'mSIRF_dot', self.handle_, ...
                other.handle_);
            sirf.Utilities.check_status('DataContainer', handle);
            re = calllib('miutilities', 'mFloatReDataFromHandle', handle);
            im = calllib('miutilities', 'mFloatImDataFromHandle', handle);
            z = complex(re, im);
            sirf.Utilities.delete(handle)
        end
        function z = uminus(self)
%***SIRF*** Overloads unary - for data containers.
            z = self*(-1);
        end
        function z = minus(self, other)
%***SIRF*** Overloads - for data containers.
%         Returns the difference of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            z = self.axpby(1, self, -1, other);
            sirf.Utilities.check_status('DataContainer:minus', z.handle_);
        end
        function z = plus(self, other)
%***SIRF*** Overloads + for data containers.
%         Returns the sum of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            z = self.axpby(1, self, 1, other);
            sirf.Utilities.check_status('DataContainer:plus', z.handle_);
        end
        function z = times(self, other)
%***SIRF*** Overloads .* for data containers.
%         Returns the elementwise product of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            z = self.same_object();
            z.handle_ = calllib('msirf', 'mSIRF_product', ...
                self.handle_, other.handle_);
            sirf.Utilities.check_status('DataContainer:times', z.handle_);
        end
        function z = rdivide(self, other)
%***SIRF*** Overloads ./ for data containers.
%         Returns the elementwise ratio of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            z = self.same_object();
            z.handle_ = calllib('msirf', 'mSIRF_ratio', ...
                self.handle_, other.handle_);
            sirf.Utilities.check_status('DataContainer:rdivide', z.handle_);
        end
        function z = mtimes(self, other)
%***SIRF*** mtimes(other) overloads * for data containers multiplication 
%         by a scalar or another data container. 
%         Returns the product self*other if other is a scalar or the dot 
%         product with other if it is a data container.
            if strcmp(class(self), class(other))
                z = self.dot(other);
                return
            end
            if isscalar(other)
                z = self.axpby(other, self, 0, self);
            else
                error('DataContainer:mtimes', ...
                    'Wrong argument type %s\n', class(other))
            end
            sirf.Utilities.check_status('DataContainer:mtimes', z.handle_);
        end
        function z = mrdivide(self, other)
%***SIRF*** mtimes(other) overloads / for data containers division
%         by a scalar. 
%         Returns the ratio self/other where other is a scalar.
            if isscalar(other)
                z = self.axpby(1.0/other, self, 0, self);
            else
                error('DataContainer:mtimes', ...
                    'Wrong argument type %s\n', class(other))
            end
            sirf.Utilities.check_status('DataContainer:mtimes', z.handle_);
        end
		function write(self, filename)
			handle = calllib('msirf', 'mSIRF_write', self.handle_, filename);
            sirf.Utilities.check_status('DataContainer:write', handle);
		end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
%***SIRF*** axpby(a, x, b, y) returns a linear combination a*x + b*y 
%         of two data containers x and y;
%         a and b: complex scalars
%         x and y: DataContainers
            sirf.Utilities.assert_validities(x, y)
            z = x.same_object();
            a = single(a);
            b = single(b);
            za = [real(a); imag(a)];
            zb = [real(b); imag(b)];
            ptr_za = libpointer('singlePtr', za);
            ptr_zb = libpointer('singlePtr', zb);
            z.handle_ = calllib('msirf', 'mSIRF_axpby', ...
                ptr_za, x.handle_, ptr_zb, y.handle_);
            sirf.Utilities.check_status('DataContainer:axpby', z.handle_);
        end
        function z = xapyb(x, a, y, b)
%***SIRF*** xapyb(x, a, y, b) returns a linear combination a*x + b*y 
%         of two data containers x and y;
%         a and b: complex scalars, or DataContainers
%         x and y: DataContainers
            sirf.Utilities.assert_validities(x, y)
            z = x.same_object();

            if isscalar(a)
                a_scalar = true;
                a = single(a);
                za = [real(a); imag(a)];
                ptr_a = libpointer('singlePtr', za);
            else
                a_scalar = false;
                sirf.Utilities.assert_validities(x, a);
                ptr_a = a.handle;
            end

            if isscalar(b)
                b_scalar = true;
                b = single(b);
                zb = [real(b); imag(b)];
                ptr_b = libpointer('singlePtr', zb);
            else
                b_scalar = false;
                sirf.Utilities.assert_validities(y, b)
                ptr_b = b.handle;
            end

            if xor(a_scalar, b_scalar)
                tmp = x.same_object();
                if a_scalar
                    tmp = b.times(y);
                    z.handle_ = calllib('msirf', 'mSIRF_axpby', ...
                        ptr_a, x.handle_, 1.0, tmp.handle_);
                else
                    tmp = a.times(x);
                    z.handle_ = calllib('msirf', 'mSIRF_axpby', ...
                        1.0, tmp.handle_, ptr_b, y.handle_);
                end
            elseif a_scalar
                z.handle_ = calllib('msirf', 'mSIRF_axpby', ...
                    ptr_a, x.handle_, ptr_b, y.handle_);
            else
                z.handle_ = calllib('msirf', 'mSIRF_xapyb', ...
                    x.handle_, ptr_a, y.handle_, ptr_b);
            end

            sirf.Utilities.check_status('DataContainer:axpby', z.handle_);
        end
    end
end
