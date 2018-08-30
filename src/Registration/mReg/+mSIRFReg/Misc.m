classdef Misc < handle
% Class for image data.

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
        function result = do_nifti_images_match(im1, im2, accuracy_percentage_of_max)
        	%Do nifti images match?
            assert(isa(im1, 'mSIRFReg.ImageData'))
            assert(isa(im2, 'mSIRFReg.ImageData'))
            h = calllib('msirfreg', 'mSIRFReg_do_nifti_images_match', im1.handle_, im2.handle_, accuracy_percentage_of_max);
            mUtilities.check_status('parameter', h)
            result = calllib('miutilities', 'mIntDataFromHandle', h);
        end
		function dump_nifti_info(to_dump)
		    %Dump metadata of one or multiple (up to 5) nifti images.
		    if ischar(to_dump)
		        h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_filename', to_dump);
		    elseif ismatrix(to_dump) && isa(to_dump, 'mSIRFReg.ImageData')
		    	if size(to_dump,2) == 1
		        	h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_im1', to_dump.handle_);
		        elseif size(to_dump,2) == 2
		        	h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_im2', to_dump(1).handle_, to_dump(2).handle_);
		        elseif size(to_dump,2) == 3
		        	h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_im3', to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_);
		        elseif size(to_dump,2) == 4
		        	h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_im4', to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_, to_dump(4).handle_);
		        elseif size(to_dump,2) == 5
		        	h = calllib('msirfreg', 'mSIRFReg_dump_nifti_info_im5', to_dump(1).handle_, to_dump(2).handle_, to_dump(3).handle_, to_dump(4).handle_, to_dump(5).handle_);
		        else
		        	error('dump_nifti_info only implemented for up to 5 images.')
		        end
		    else
		    	error('dump_nifti_info requires filename, SIRFImageData or a list of SIRFImageData.')
		    end
		    mUtilities.check_status('parameter', h)
		end
    end
end