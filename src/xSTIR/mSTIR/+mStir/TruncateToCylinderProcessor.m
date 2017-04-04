classdef TruncateToCylinderProcessor < mStir.ImageDataProcessor
% Class for the image filter that zeroes the image outside the cylinder
% of the same xy-diameter and z-size as those of the image.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
        function self = TruncateToCylinderProcessor()
            self.name = 'TruncateToCylindricalFOVImageProcessor';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_strictly_less_than_radius(self, flag)
%***SIRF*** set_strictly_less_than_radius(flag) specifies whether the area
%         not affected by filtering is strictly inside the cylinder 
%         (flag = True) or not (flag = False).
            if flag
                str = 'true';
            else
                str = 'false';
            end
            mStir.setParameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', str, 'c')
        end
        function flag = get_strictly_less_than_radius(self)
%***SIRF*** Returns the answer to the question: Is the area not affected by 
%         filtering strictly inside the cylinder?
            flag = mStir.parameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', 'i');
        end
    end
end