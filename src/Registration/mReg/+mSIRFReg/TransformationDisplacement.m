classdef TransformationDisplacement < mSIRFReg.Transformation
% Class for displacement transformations.

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
            name = 'SIRFRegTransformationDisplacement';
        end
    end
    methods
        function self = TransformationDisplacement(src)
            narginchk(0,1)
            self.name = 'SIRFRegTransformationDisplacement';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, src);
            elseif isa(src, 'mSIRFReg.NiftiImage3DDisplacement')
                self.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegTransformationDisplacement_construct_from_NiftiImage3DDisplacement', src.handle_);
            else
                error('TransformationDisplacement accepts no args, filename or NiftiImage3DDisplacement.')
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
    end
end