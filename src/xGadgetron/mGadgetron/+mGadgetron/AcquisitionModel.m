classdef AcquisitionModel < handle
% Class for MR acquisition model, an operator that maps images into
% simulated acquisitions.

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

    properties
        handle_
        name_
    end
    methods
        function self = AcquisitionModel(acqs, imgs)
            self.name_ = 'MR_AcquisitionModel';
            self.handle_ = calllib('mgadgetron', 'mGT_AcquisitionModel',...
                acqs.handle_, imgs.handle_);
            mUtil.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function set_coil_sensitivity_maps(self, csms)
%         Specifies the coil sensitivity maps to be used by the model.
%         csm: CoilSensitivityData
            handle = calllib('mgadgetron', 'mGT_setCSMs', ...
                self.handle_, csms.handle_);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function acqs = forward(self, imgs)
%         Projects an image into (simulated) acquisitions space.
%         The resulting acquisition data simulates the actual data
%         expected to be received from the scanner.
%         image: ImageData
            acqs = mGadgetron.AcquisitionData();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelForward', ...
                self.handle_, imgs.handle_);
            mUtil.checkExecutionStatus(self.name_, acqs.handle_);
        end
        function imgs = backward(self, acqs)
%         Back-projects acquisition data into image space using a complex
%         transpose of the forward projection.
%         ad: AcquisitionData
            imgs = mGadgetron.ImageData();
            imgs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelBackward', ...
                self.handle_, acqs.handle_);
            mUtil.checkExecutionStatus(self.name_, imgs.handle_);
        end
    end
end