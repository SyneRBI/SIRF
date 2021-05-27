classdef ImageData < sirf.SIRF.DataContainer
% INTERNAL USE ONLY.
% Class for an abstract data container.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
% Copyright 2018 - 2020 University College London
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

	methods
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function same = eq(self, other)
            assert(isa(other, 'sirf.SIRF.ImageData'));
            h = calllib('msirf', 'mSIRF_equalImages', self.handle_, other.handle_);
            sirf.Utilities.check_status('ImageData:eq', h);
            same = logical(calllib('miutilities', 'mIntDataFromHandle', h));
            sirf.Utilities.delete(h)
        end
        function diff = ne(self, other)
        	diff = ~(self == other);
        end
        function read(self, filename, engine, verb)
            self.handle_ = calllib('msirf', 'mSIRF_readImageData', filename, engine, verb);
            sirf.Utilities.check_status('ImageData:read', self.handle_);
        end
        function fill(self, image)
            calllib('msirf', 'mSIRF_fillImageFromImage', self.handle_, image.handle_);
        end
		function geom_info = get_geometrical_info(self)
			% Get the image's geometrical info.
			geom_info = sirf.SIRF.GeometricalInfo();
			geom_info.handle_ = calllib('msirf', 'mSIRF_ImageData_get_geom_info', self.handle_);
		end
		function reorient(self, geom_info)
            % Reorient image. Requires that dimensions match.
            assert(isa(geom_info, 'sirf.SIRF.GeometricalInfo'));
            h = calllib('msirf', 'mSIRF_ImageData_reorient', self.handle_, geom_info.handle_);
            sirf.Utilities.check_status([self.name ':reorient'], h);
            sirf.Utilities.delete(h)
        end
	end
end
