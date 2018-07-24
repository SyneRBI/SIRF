classdef AcquisitionModel < handle
%     Class for PET acquisition model objects.
%     Relates an image x to the simulated acquisition data y as
%     (F)    y = [1/n](G x + [a]) + [b]
%     where:
%     G is the geometric (ray tracing) projector from the image voxels
%     to the scanner's pairs of detectors (bins);
%     a and b are optional additive and background terms representing
%     the effects of accidental coincidences and scattering; assumed to be 
%     0 if not present;
%     n is an optional bin normalization term representing the inverse of
%     detector (bin) efficiencies; assumed to be 1 if not present.
%     The computation of y for a given x by the above formula (F) is
%     referred to as (forward) projection, and the computation of
%     (B)    z = G' [1/n] y
%     where G' is the transpose of G, is referred to as backprojection.

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

    properties
        name
        handle_
    end
    methods (Static)
        function name = class_name()
            name = 'AcquisitionModel';
        end
    end
    methods
        function self = AcquisitionModel()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_additive_term(self, at)
%***SIRF*** sets the additive term a in (F);
%         at:  an AcquisitionData object containing a.
            mUtilities.assert_validity(at, 'AcquisitionData')
            mSTIR.setParameter(self.handle_, 'AcquisitionModel', ...
                'additive_term', at, 'h');
        end
        function set_background_term(self, bt)
%***SIRF*** sets the background term b in (F);
%         bt:  an AcquisitionData object containing b.
            mUtilities.assert_validity(bt, 'AcquisitionData')
            mSTIR.setParameter(self.handle_, 'AcquisitionModel', ...
                'background_term', bt, 'h');
        end
        function set_acquisition_sensitivity(self, asm)
%***SIRF*** sets the acquisition sensitivity model responsible for n in (F);
%         asm: an AcquisitionSensitivityModel object.
            mUtilities.assert_validity(asm, 'AcquisitionSensitivityModel')
            mSTIR.setParameter(self.handle_, 'AcquisitionModel', ...
                'asm', asm, 'h');
        end
        function set_up(self, acq_templ, img_templ)
%***SIRF*** sets up the object with appropriate geometric information.
%         This function needs to be called before performing forward- or 
%         backprojections;
%         Usage: 
%             set_up(acq_templ, img_templ);
%         acq_templ:  an AcquisitionData object used as a template for
%                     creating an AcquisitionData object to store forward
%                     projection;
%         img_templ:  an ImageData object used as a template for creating an
%                     ImageData object to store backprojection.
            mUtilities.assert_validity(acq_templ, 'AcquisitionData')
            mUtilities.assert_validity(img_templ, 'ImageData')
            h = calllib...
                ('mstir', 'mSTIR_setupAcquisitionModel',...
                self.handle_, acq_templ.handle_, img_templ.handle_);
            mUtilities.check_status([self.name ':set_up'], h)
            mUtilities.delete(h)
            %calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function ad = forward(self, image, subset_num, num_subsets)
%***SIRF*** computes the forward projection of ImageData x = image given by 
%         (F) above in the main class documentation.
%         Usage: 
%             acq_data = forward(image);
%         image:  an ImageData object containing x;
            mUtilities.assert_validity(image, 'ImageData')
            if nargin < 4
                subset_num = 0;
                num_subsets = 1;
            end
            ad = mSTIR.AcquisitionData();
            ad.handle_ = calllib('mstir', 'mSTIR_acquisitionModelFwd',...
                self.handle_, image.handle_, subset_num, num_subsets);
            mUtilities.check_status([self.name ':forward'], ad.handle_)
        end
        function image = backward(self, ad, subset_num, num_subsets)
%***SIRF*** returns the backprojection z for y = ad given by (B) above;
%         ad:  an AcquisitionData object containing y.
            mUtilities.assert_validity(ad, 'AcquisitionData')
            if nargin < 4
                subset_num = 0;
                num_subsets = 1;
            end
            image = mSTIR.ImageData();
            image.handle_ = calllib('mstir', 'mSTIR_acquisitionModelBwd',...
                self.handle_, ad.handle_, subset_num, num_subsets);
            mUtilities.check_status...
                ([self.name ':backward'], image.handle_)
        end
    end
end