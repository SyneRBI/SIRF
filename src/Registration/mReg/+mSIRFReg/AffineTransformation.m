classdef AffineTransformation < mSIRFReg.Transformation
% Class for affine transformations.

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
            name = 'SIRFRegAffineTransformation';
        end
    end
    methods
        function self = AffineTransformation(src)
            narginchk(0,1)
            self.name = 'SIRFRegAffineTransformation';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, src);
            elseif isnumeric(src)
                assert(all(size(src)==[4, 4]))
                ptr_v = libpointer('singlePtr', single(src));
                self.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_construct_from_TM', ptr_v);
            else
                error('AffineTransformation accepts no args, filename or 4x4 array.')
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function value = eq(self, other)
            %Overload comparison operator.
            assert(isa(other, 'mSIRFReg.AffineTransformation'))
            handle = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_equal', self.handle_, other.handle_);
            mUtilities.check_status('SIRFRegAffineTransformation:eq', handle);
            value = logical(calllib('miutilities', 'mIntDataFromHandle', handle));
            mUtilities.delete(handle)
        end
        function value = ne(self, other)
            %Overload comparison operator.
            value = (self == other);
        end
        function mat = mtimes(self, other)
            %Overload multiplication operator.
            assert(isa(other, 'mSIRFReg.AffineTransformation'))
            mat = mSIRFReg.AffineTransformation();
            mat.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_equal', self.handle_, other.handle_);
            mUtilities.check_status('SIRFRegAffineTransformation:mtimes', mat);
        end

        function mat = deep_copy(self)
            %Deep copy.
            mat = mSIRFReg.AffineTransformation();
            mat.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_deep_copy', self.handle_);
            mUtilities.check_status('SIRFRegAffineTransformation:mtimes', mat);
        end

        function save_to_file(self, filename)
            %Save to file.
            calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_save_to_file', self.handle_, filename);
        end
        function value = get_determinant(self)
            %Get determinant.
            value = mSTIR.parameter(self.handle_, self.name, 'determinant', 'f');
        end
        function tm = as_array(self)
            %Get forward transformation matrix.
            ptr_v = libpointer('singlePtr', zeros(4, 4));
            calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_as_array', self.handle_, ptr_v);
            tm = ptr_v.Value;
        end    
        function tm = get_inverse(self)
            %Get forward transformation matrix.
            tm = mSIRFReg.AffineTransformation();
            tm.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_get_inverse', self.handle_);
            mUtilities.check_status('SIRFRegAffineTransformation:get_inverse', tm);
        end
    end
    methods(Static)
        function mat = get_identity()
            %Get identity matrix.
            mat = mSIRFReg.AffineTransformation();
            mat.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegAffineTransformation_get_identity');
        end
    end
end
