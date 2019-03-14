classdef TruncateToCylinderProcessor < sirf.STIR.ImageDataProcessor
% Class for the image filter that zeroes the image outside a cylinder.

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
        function self = TruncateToCylinderProcessor()
%         Creates a TruncateToCylinderProcessor object.
%         The application of this processor to an image zeroes its values
%         outside the vertical cylinder inscribed into the image's bounding 
%         box. The treatment of values on the cylinder boundary is
%         defined by set_strictly_less_than_radius() method.
            self.name_ = 'TruncateToCylindricalFOVImageProcessor';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_strictly_less_than_radius(self, flag)
%***SIRF*** Defines the filter behaviour on the boundary of the cylinder.
%         set_strictly_less_than_radius(flag) specifies whether the area
%         not affected by filtering is strictly inside the cylinder 
%         (flag = True) or not (flag = False).
            if flag
                str = 'true';
            else
                str = 'false';
            end
            sirf.STIR.setParameter(self.handle_,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', str, 'c')
        end
        function flag = get_strictly_less_than_radius(self)
%***SIRF*** Exposes the filter behaviour on the boundary of the cylinder.
%         Returns the answer to the question: Is the area not affected by 
%         filtering strictly inside the cylinder?
            flag = sirf.STIR.parameter(self.handle_,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', 'i');
        end
    end
end