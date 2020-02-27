classdef SPMRegistration < sirf.Reg.NiftyRegistration
% Registration class using SPM.

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

    methods(Static)
        function name = class_name()
            name = 'SPMRegistration';
        end
    end
    methods
        function self = SPMRegistration()
            self.name = 'SPMRegistration';
            self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            sirf.Utilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function tm = get_transformation_matrix_forward(self, idx)
            %Get forward transformation matrix. 1-based.
            if nargin < 2; idx=1; end
            tm = sirf.Reg.AffineTransformation();
            tm.handle_ = calllib('mreg', 'mReg_SPMRegistration_get_TM', self.handle_, 'forward', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_transformation_matrix_forward'], tm.handle_);
        end
        function tm = get_transformation_matrix_inverse(self, idx)
            %Get inverse transformation matrix. 1-based.
            if nargin < 2; idx=1; end
            tm = sirf.Reg.AffineTransformation();
            tm.handle_ = calllib('mreg', 'mReg_SPMRegistration_get_TM', self.handle_, 'inverse', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_transformation_matrix_inverse'], tm.handle_);
        end
        function set_working_folder(self, working_folder)
            %Set working folder.
            sirf.Reg.setParameter(self.handle_, 'SPMRegistration', 'working_folder', working_folder, [])
        end
        function set_working_folder_file_overwrite(self, working_folder_file_overwrite)
            %Set file overwrite in working folder.
            if working_folder_file_overwrite
                working_folder_file_overwrite = 1;
            else
                working_folder_file_overwrite = 0;
            end
            sirf.Reg.setParameter(self.handle_, 'SPMRegistration', 'working_folder_file_overwrite', working_folder_file_overwrite, 'i')
        end
        function set_delete_temp_files(self, delete_temp_files)
            %Set delete temp files.
            if delete_temp_files
                delete_temp_files = 1;
            else
                delete_temp_files = 0;
            end
            sirf.Reg.setParameter(self.handle_, 'SPMRegistration', 'delete_temp_files', delete_temp_files, 'i')
        end
    end
end
