classdef NiftiImageData3DDeformation < mReg.NiftiImageData3DTensor & mReg.Transformation
% Class for deformation image data.

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

    methods(Static)
        function name = class_name()
            name = 'NiftiImageData3DDeformation';
        end
        function obj = same_object()
            obj = mReg.NiftiImageData3DDeformation();
        end
    end
    methods
        function self = NiftiImageData3DDeformation(src1, src2, src3)
            narginchk(0,3)
            self.name = 'NiftiImageData3DDeformation';
            if nargin < 1
                self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            elseif ischar(src1)
                self.handle_ = calllib('mreg', 'mReg_objectFromFile', self.name, src1);
            elseif nargin == 3 && isa(src1, 'mSIRF.ImageData') && isa(src2, 'mSIRF.ImageData') && isa(src3, 'mSIRF.ImageData')
                self.handle_ = calllib('mreg', 'mReg_NiftiImageData3DTensor_construct_from_3_components', self.name, src1.handle_, src2.handle_, src3.handle_);                
            elseif isa(src1, 'mReg.NiftiImageData3DDisplacement')
                self.handle_ = calllib('mreg', 'mReg_NiftiImageData3DDeformation_create_from_disp', src1.handle_);
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
    methods(Static)
        function z = compose_single_deformation(trans, ref)
	    	%Compose up to transformations into single deformation.
		    assert(isa(ref, 'mSIRF.ImageData'))
		    assert(isa(trans, 'mReg.Transformation'))
		    if isrow(trans)
                trans=trans'; 
            end
		    assert(iscolumn(trans));
            num_trans = size(trans,1);
		    if num_trans == 1
		    	z = trans(1).get_as_deformation_field(ref);
		        return
            end
            % This is ugly. Store each type in a single string (need to do this because I can't get
            % virtual methods to work for multiple inheritance (deformation/displacement are both
            % nifti images and transformations).
            types = '';
            for n = 1:num_trans
                if isa(trans(n),'mReg.AffineTransformation')
                    types = [types '1'];
                elseif isa(trans(n),'mReg.NiftiImageData3DDisplacement')
                    types = [types '2'];
                elseif isa(trans(n),'mReg.NiftiImageData3DDeformation')
                    types = [types '3'];
                end
            end
		    z = mReg.NiftiImageData3DDeformation();
		    if num_trans == 2
		        z.handle_ = calllib('mreg', 'mReg_NiftiImageData3DDeformation_compose_single_deformation',...
		        	ref.handle_, num_trans, types, trans(1).handle_, trans(2).handle_, [], [], []);
		    elseif num_trans == 3
		        z.handle_ = calllib('mreg', 'mReg_NiftiImageData3DDeformation_compose_single_deformation',...
		            ref.handle_, num_trans, types, trans(1).handle_, trans(2).handle_, trans(3).handle_, [], []);
		    elseif num_trans == 4
		        z.handle_ = calllib('mreg', 'mReg_NiftiImageData3DDeformation_compose_single_deformation',...
		            ref.handle_, num_trans, types, trans(1).handle_, trans(2).handle_, trans(3).handle_, trans(4).handle_, []);
		    elseif num_trans == 5
		        z.handle_ = calllib('mreg', 'mReg_NiftiImageData3DDeformation_compose_single_deformation',...
		            ref.handle_, num_trans, types, trans(1).handle_, trans(2).handle_, trans(3).handle_, trans(4).handle_, trans(5).handle_);
		    else
		        error('compose_transformations_into_single_deformation only implemented for up to 5 transformations.')
		    end
		    mUtilities.check_status('compose_transformations_into_single_deformation', z.handle_);
		end
    end
end