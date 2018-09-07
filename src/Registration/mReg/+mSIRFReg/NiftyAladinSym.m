classdef NiftyAladinSym < mSIRFReg.SIRFReg
% Registration class using NiftyReg's symmetric aladin.

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

    methods(Static)
        function name = class_name()
            name = 'SIRFRegNiftyAladinSym';
        end
    end
    methods
        function self = NiftyAladinSym()
            self.name = 'SIRFRegNiftyAladinSym';
            self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function save_transformation_matrix_fwrd(self, filename)
            %Save forward transformation matrix.
            narginchk(2,2);
            assert(ischar(filename))
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix', self.handle_, filename, 'fwrd');
            mUtilities.check_status([self.name ':save_transformation_matrix_fwrd'], h);
            mUtilities.delete(h)
        end
        function save_transformation_matrix_back(self, filename)
            %Save backward transformation matrix.
            narginchk(2,2);
            assert(ischar(filename))
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix', self.handle_, filename, 'back');
            mUtilities.check_status([self.name ':save_transformation_matrix_back'], h);
            mUtilities.delete(h)
        end
        function tm = get_transformation_matrix_fwrd(self)
            %Get forward transformation matrix.
            ptr_v = libpointer('singlePtr', zeros(4, 4));
            calllib('msirfreg', 'mSIRFReg_SIRFReg_get_TM', self.handle_, ptr_v, 'fwrd');
            tm = ptr_v.Value;
        end
        function tm = get_transformation_matrix_back(self)
            %Get backwards transformation matrix.
            ptr_v = libpointer('singlePtr', zeros(4, 4));
            calllib('msirfreg', 'mSIRFReg_SIRFReg_get_TM', self.handle_, ptr_v, 'back');
            tm = ptr_v.Value;
        end
    end
end