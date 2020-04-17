classdef GeometricalInfo < handle
% INTERNAL USE ONLY.
% Class for image geometrical info.

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
            name = 'GeometricalInfo';
        end
        function obj = same_object()
            obj = sirf.Reg.GeometricalInfo();
        end
    end
    methods
        function self = GeometricalInfo()
        	self.handle_ = [];
        end
        function print_info(self)
        	% Print geometrical information
        	h = calllib('msirf', 'mSIRF_GeomInfo_print', self.handle_);
            sirf.Utilities.check_status([self.name ':print_info'], h);
            sirf.Utilities.delete(h)
        end
        function value = get_offset(self)
        	% Offset is the LPS coordinate of the centre of the first voxel.
        	ptr_i = libpointer('int32Ptr', zeros(1, 3));
            calllib('msirf', 'mSIRF_GeomInfo_get_offset', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function value = get_spacing(self)
        	% Spacing is the physical distance between voxels in each dimension.
        	ptr_i = libpointer('singlePtr', zeros(1, 3));
            calllib('msirf', 'mSIRF_GeomInfo_get_spacing', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function value = get_size(self)
        	% Size is the number of voxels in each dimension.
        	ptr_i = libpointer('int32Ptr', zeros(1, 3));
            calllib('msirf', 'mSIRF_GeomInfo_get_size', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function value = get_direction_matrix(self)
        	% Each row gives a vector dictating the direction of the axis in LPS physical space.
        	ptr_i = libpointer('singlePtr', zeros(3, 3));
            calllib('msirf', 'mSIRF_GeomInfo_get_direction_matrix', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
        function value = get_index_to_physical_point_matrix(self)
        	% Get the 4x4 affine matrix that converts an index to a point in LPS physical space.
        	ptr_i = libpointer('singlePtr', zeros(4,4));
            calllib('msirf', 'get_index_to_physical_point_matrix', self.handle_, ptr_i);
            value = ptr_i.Value;
        end
    end
end