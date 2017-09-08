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
        same_object(self)
    end
    methods
        function self = DataContainer()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
%         function num = number(self)
% %***SIRF*** Returns the number of items in this container.
%             handle = calllib('mgadgetron', 'mGT_dataItems', self.handle_);
%             mUtilities.check_status('DataContainer', handle);
%             num = calllib('miutilities', 'mIntDataFromHandle', handle);
%             mUtilities.delete(handle)
%         end
        function r = norm(self)
%***SIRF*** Returns the 2-norm of this data container viewed as a vector.
            handle = calllib('mstir', 'mSTIR_norm', self.handle_);
            mUtilities.check_status('DataContainer', handle);
            r = calllib('miutilities', 'mFloatDataFromHandle', handle);
            mUtilities.delete(handle)
        end
        function r = dot(self, other)
%***SIRF*** Returns the dot product of this data container with another one 
%         viewed as vectors.
            handle = calllib('mstir', 'mSTIR_dot', self.handle_, ...
                other.handle_);
            mUtilities.check_status('DataContainer', handle);
            r = calllib('miutilities', 'mFloatDataFromHandle', handle);
            mUtilities.delete(handle)
        end
        function z = plus(self, other)
%***SIRF*** Overloads - for data containers.
%         Returns the difference of this data container with another one
%         viewed as vectors.
            z = self.same_object();
            z.handle_ = calllib('mstir', 'mSTIR_axpby', ...
                1.0, self.handle_, 1.0, other.handle_);
        end
        function z = minus(self, other)
%***SIRF*** Overloads - for data containers.
%         Returns the difference of this data container with another one
%         viewed as vectors.
            z = self.same_object();
            z.handle_ = calllib('mstir', 'mSTIR_axpby', ...
                1.0, self.handle_, -1.0, other.handle_);
        end
        function z = mtimes(self, other)
%***SIRF*** mtimes(other) overloads * for data containers multiplication 
%         by a scalar or another data container. 
%         Returns the product self*other if other is a scalar or the dot 
%         product with other if it is a data container.
            if isobject(other)
                z = self.dot(other);
            elseif isreal(other)
                z = self.same_object();
                z.handle_ = calllib('mstir', 'mSTIR_mult', other, self.handle_);
            else
                error('DataContainer:mtimes', 'Wrong multiplier');
            end
        end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
%***SIRF*** axpby(a, x, b, y) returns a linear combination a*x + b*y 
%         of two data containers x and y;
%         a and b: real scalars
%         x and y: DataContainers
            z = self.same_object();
            z.handle_ = calllib('mstir', 'mSTIR_axpby', ...
                a, x.handle_, b, y.handle_);
        end
    end
end