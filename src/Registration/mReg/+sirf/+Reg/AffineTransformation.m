classdef AffineTransformation < sirf.Reg.Transformation
% Class for affine transformations.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018-2020 University College London
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
            name = 'AffineTransformation';
        end
    end
    methods
        function self = AffineTransformation(src1,src2)
            narginchk(0,2)
            self.name = 'AffineTransformation';
            if nargin < 1
                self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            elseif ischar(src1)
                self.handle_ = calllib('mreg', 'mReg_objectFromFile', self.name, src1);
            elseif isnumeric(src1) && nargin==1
                assert(all(size(src1)==[4, 4]))
                ptr_v = libpointer('singlePtr', single(src1));
                self.handle_ = calllib('mreg', 'mReg_AffineTransformation_construct_from_TM', ptr_v);
            elseif isnumeric(src1) && isa(src2, 'sirf.Reg.Quaternion')
                ptr_v = libpointer('singlePtr', single(src1));
                self.handle_ = calllib('mreg', 'mReg_AffineTransformation_construct_from_trans_and_quaternion', ptr_v, src2.handle_);
            elseif isnumeric(src1) && isnumeric(src2)
                ptr1_v = libpointer('singlePtr', single(src1));
                ptr2_v = libpointer('singlePtr', single(src2));
                self.handle_ = calllib('mreg', 'mReg_AffineTransformation_construct_from_trans_and_euler', ptr1_v, ptr2_v);
            else
                error('AffineTransformation accepts no args, filename, 4x4 array or translation with quaternion.')
            end
            sirf.Utilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function value = eq(self, other)
            %Overload comparison operator.
            assert(isa(other, 'sirf.Reg.AffineTransformation'))
            handle = calllib('mreg', 'mReg_AffineTransformation_equal', self.handle_, other.handle_);
            sirf.Utilities.check_status('AffineTransformation:eq', handle);
            value = logical(calllib('miutilities', 'mIntDataFromHandle', handle));
            sirf.Utilities.delete(handle)
        end
        function value = ne(self, other)
            %Overload comparison operator.
            value = ~(self == other);
        end
        function mat = mtimes(self, other)
            %Overload multiplication operator.
            assert(isa(other, 'sirf.Reg.AffineTransformation'))
            mat = sirf.Reg.AffineTransformation();
            mat.handle_ = calllib('mreg', 'mReg_AffineTransformation_mul', self.handle_, other.handle_);
            sirf.Utilities.check_status('AffineTransformation:mtimes', mat.handle_);
        end

        function mat = deep_copy(self)
            %Deep copy.
            mat = sirf.Reg.AffineTransformation();
            mat.handle_ = calllib('mreg', 'mReg_AffineTransformation_deep_copy', self.handle_);
            sirf.Utilities.check_status('AffineTransformation:mtimes', mat.handle_);
        end

        function write(self, filename)
            %Save to file.
            calllib('mreg', 'mReg_AffineTransformation_write', self.handle_, filename);
        end
        function value = get_determinant(self)
            %Get determinant.
            value = sirf.Reg.parameter(self.handle_, self.name, 'determinant', 'f');
        end
        function tm = as_array(self)
            %Get forward transformation matrix.
            ptr_v = libpointer('singlePtr', zeros(4, 4));
            calllib('mreg', 'mReg_AffineTransformation_as_array', self.handle_, ptr_v);
            tm = ptr_v.Value';
        end    
        function tm = get_inverse(self)
            %Get forward transformation matrix.
            tm = sirf.Reg.AffineTransformation();
            tm.handle_ = calllib('mreg', 'mReg_AffineTransformation_get_inverse', self.handle_);
            sirf.Utilities.check_status('AffineTransformation:get_inverse', tm.handle_);
        end
        function eul = get_Euler_angles(self)
            %Get Euler angles of transformation matrix (XYZ).
            ptr_v = libpointer('singlePtr', zeros(1,3));
            calllib('mreg', 'mReg_AffineTransformation_get_Euler_angles', self.handle_, ptr_v);
            eul = ptr_v.Value;
        end
        function quat = get_quaternion(self)
            %Get quaternion.
            quat = sirf.Reg.Quaternion([0,0,0,0]);
            quat.handle_ = calllib('mreg', 'mReg_AffineTransformation_get_quaternion', self.handle_);
            sirf.Utilities.check_status('AffineTransformation:get_quaternion', quat.handle_);
        end
    end
    methods(Static)
        function mat = get_identity()
            %Get identity matrix.
            mat = sirf.Reg.AffineTransformation();
            mat.handle_ = calllib('mreg', 'mReg_AffineTransformation_get_identity');
        end
        function tm = get_average(to_average)
            %Average transformation matrices.
            vec = sirf.SIRF.DataHandleVector();
            for n = 1:length(to_average)
                vec.push_back(to_average(n).handle_);
            end
            tm = sirf.Reg.AffineTransformation();
            tm.handle_ = calllib('mreg', 'mReg_AffineTransformation_get_average', vec.handle_);
            sirf.Utilities.check_status('AffineTransformation:get_average', tm.handle_)
        end
    end
end
