classdef ImageData < sirf.SIRF.DataContainer
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

	methods
		function geom_info = get_geometrical_info(self)
			% Get the image's geometrical info.
			geom_info = sirf.SIRF.GeometricalInfo();
			geom_info.handle_ = calllib('msirf', 'mSIRF_ImageData_get_geom_info', self.handle_);
		end
		function reorient(self, geom_info)
            % Reorient image. Requires that dimensions and spacing match.
            assert(isa(geom_info, 'sirf.SIRF.GeometricalInfo'));
            h = calllib('msirf', 'mSIRF_ImageData_reorient', self.handle_, geom_info.handle_);
            sirf.Utilities.check_status([self.name ':reorient'], h);
            sirf.Utilities.delete(h)
        end
	end
end