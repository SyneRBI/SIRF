classdef DataContainer < handle
% INTERNAL USE ONLY.
% Class for an abstract data container.

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
        function out = times(self, other, out)
%***SIRF*** Overloads .* for data containers.
%         Returns the elementwise product of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            if nargin==2
                out = self.copy();
            end
            h = calllib('msirf', 'mSIRF_multiply', ...
                self.handle_, other.handle_, out.handle_);
            sirf.Utilities.check_status('DataContainer:times', h);
            sirf.Utilities.delete(h)
        end
        function z = rdivide(self, other)
%***SIRF*** Overloads ./ for data containers.
%         Returns the elementwise ratio of this data container with another one
%         viewed as vectors.
            sirf.Utilities.assert_validities(self, other)
            if nargin==2
                out = self.copy();
            end
            h = calllib('msirf', 'mSIRF_divide', ...
                self.handle_, other.handle_, out.handle_);
            sirf.Utilities.check_status('DataContainer:divide', h);
            sirf.Utilities.delete(h)
        end
        function z = mtimes(self, other)
%***SIRF*** mtimes(other) overloads * for data containers multiplication 
%         by a scalar or another data container. 
%         Returns the product self*other if other is a scalar or the dot 
%         product with other if it is a data container.
            %if isobject(other)
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
%***SIRF*** mtimes(other) overloads / for data containers multiplication 
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
            %assert(strcmp(class(x), class(y)))
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
    end
end