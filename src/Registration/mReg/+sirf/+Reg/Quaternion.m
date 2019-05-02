classdef Quaternion < handle
% Class for quaternions.

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
            name = 'Quaternion';
        end
    end
    methods
        function self = Quaternion(src)
            narginchk(1,1)
            self.name = 'Quaternion';
            if isnumeric(src)
                assert(length(src)==4,'Quaternion should consist of 4 values');
                ptr_v = libpointer('singlePtr', single(src));
                self.handle_ = calllib('mreg', 'mReg_Quaternion_construct_from_array', ptr_v);
            elseif isa(src, 'sirf.Reg.AffineTransformation')
                self.handle_ = calllib('mreg', 'mReg_Quaternion_construct_from_AffineTransformation', src.handle_);
            else
                error('AffineTransformation accepts no args, filename or array.')
            end
            sirf.Utilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function arr = as_array(self)
            %Get quaternion as array.
            ptr_v = libpointer('singlePtr', zeros(1, 4));
            calllib('mreg', 'mReg_Quaternion_as_array', self.handle_, ptr_v);
            arr = ptr_v.Value;
        end
    end
    methods(Static)
        function quat = get_average(to_average)
            %Get average of quaternions.
            assert(ismatrix(to_average) && isa(to_average, 'sirf.Reg.Quaternion'), 'Quaternion.get_average: give list of Quaternion.')
            vec = sirf.SIRF.DataHandleVector();
            for n = 1:length(to_average)
                vec.push_back(to_average(n).handle_);
            end
            quat = sirf.Reg.Quaternion([0,0,0,0]);
            quat.handle_ = calllib('mreg', 'mReg_Quaternion_get_average', vec.handle_);
            sirf.Utilities.check_status('Quaternion:get_average', quat.handle_)
        end
    end
end
