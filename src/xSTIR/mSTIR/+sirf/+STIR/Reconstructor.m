classdef Reconstructor < handle
% Class for a generic PET reconstructor object.

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

    properties (Constant)
        R = 'Reconstruction';
    end
    properties
        handle_
        input
        image
    end
    methods
        function self = Reconstructor()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mstir', 'mDeleteDataHandle', self.handle_)
            end
        end
        function set_input(self, input_data)
            mUtilities.assert_validity(input_data, 'AcquisitionData')
            sirf.STIR.setParameter...
                (self.handle_, self.R, 'input_data', input_data, 'h')
        end
        function process(self)
%***SIRF*** Reconstruct the image 
%         by applying currently set range of
%         iterations to the current image estimate.
            if isempty(self.image)
                error('Reconstructor:process', 'current estimate not set')
            end
            mUtilities.assert_validity(self.image, 'ImageData')
            h = calllib('mstir', 'mSTIR_runReconstruction',...
                self.handle_, self.image.handle_);
            mUtilities.check_status('Reconstructor:process', h)
            mUtilities.delete(h)
        end
        function image = get_output(self)
            image = self.image;
        end
        function reconstruct(self, image)
%***SIRF*** Reconstruct the image 
%         by applying currently set range of
%         iterations to the image estimate specified by the argument.
            mUtilities.assert_validity(image, 'ImageData')
            h = calllib('mstir', 'mSTIR_runReconstruction',...
                self.handle_, image.handle_);
            mUtilities.check_status([self.IR ':reconstruct'], h)
            mUtilities.delete(h)
        end
        function set_output_filename_prefix(self, prefix)
%***SIRF*** Specifies the naming for the output files.
%         This method sets the prefix that is used to determine the
%         filename for output files with the image estimates at
%         different sub-iterations.
%         Usage: 
%             recon.set_output_filename_prefix(prefix);
%         prefix: string with prefix of filename.
%                 Each file will be named [prefix '_' subiter_num], 
%                 where subiter_num is the number of the sub-iteration 
%                 at which the respective image estimate was saved.
            sirf.STIR.setParameter(self.handle_, self.R, 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end