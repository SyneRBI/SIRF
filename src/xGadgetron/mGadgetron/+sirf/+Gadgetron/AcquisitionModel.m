classdef AcquisitionModel < handle
% Class for MR acquisition model.
% The MR acquisition model describes the transformation (encoding) 
% of image data into MR acquisitions 

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
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

    properties
        handle_
        name_
    end
    methods
        function self = AcquisitionModel(acq_template, img_template)
 %        AcquisitionModel(acq_templ, img_templ) creates an MR acquisition model 
 %        based on two specified templates of types AcquisitionData and ImageData 
 %        respectively.
            self.name_ = 'MR_AcquisitionModel';
            sirf.Utilities.assert_validity(acq_template, 'AcquisitionData')
            sirf.Utilities.assert_validity(img_template, 'ImageData')
            self.handle_ = calllib('mgadgetron', 'mGT_AcquisitionModel',...
                acq_template.handle_, img_template.handle_);
            sirf.Utilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function set_coil_sensitivity_maps(self, csms)
%***SIRF*** Instructs to use the coil sensitivity maps specified by the argument
%         (a CoilSensitivityData object).
            sirf.Utilities.assert_validity(csms, 'CoilSensitivityData')
            handle = calllib('mgadgetron', 'mGT_setCSMs', ...
                self.handle_, csms.handle_);
            sirf.Utilities.check_status(self.name_, handle);
            sirf.Utilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function acqs = forward(self, image)
%***SIRF*** Returns the forward projection of the specified ImageData argument
%         simulating the actual data expected to be received from the scanner.
            sirf.Utilities.assert_validity(image, 'ImageData')
            acqs = sirf.Gadgetron.AcquisitionData();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelForward', ...
                self.handle_, image.handle_);
            sirf.Utilities.check_status(self.name_, acqs.handle_);
        end
        function imgs = backward(self, acqs)
%***SIRF*** Backprojects the acquisition data specified by the argument
%         of AcquisitionData type into image space using a complex
%         transpose of the forward projection.
            sirf.Utilities.assert_validity(acqs, 'AcquisitionData')
            imgs = sirf.Gadgetron.ImageData();
            imgs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelBackward', ...
                self.handle_, acqs.handle_);
            sirf.Utilities.check_status(self.name_, imgs.handle_);
        end
    end
end
