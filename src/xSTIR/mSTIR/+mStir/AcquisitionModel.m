classdef AcquisitionModel < handle
%     Class for a PET acquisition model that relates an image x to the
%     acquisition data y as
%     (F)    y = [1/n](G x + [a]) + [b]
%     where:
%     G is the geometric (ray tracing) projector from the image voxels
%     to the scanner's pairs of detectors (bins);
%     a and b are otional additive and background terms representing
%     the effects of noise and scattering; assumed to be 0 if not present;
%     n is an optional bin normalization term representing the inverse of
%     detector (bin) efficiencies; assumed to be 1 if not present.
%     The computation of y for a given x by the above formula (F) is
%     referred to as forward projection, and the computation of
%     (B)    z = G' m y
%     where G' is the transpose of G and m = 1/n, is referred to as 
%     backward projection.

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
        name
        handle
    end
    methods
        function self = AcquisitionModel()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_additive_term(self, at)
%***SIRF*** Sets the additive term a in (F);
%         at:  an AcquisitionData object containing a.
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', at, 'h');
        end
        function set_background_term(self, bt)
%***SIRF*** Sets the background term b in (F);
%         bt:  an AcquisitionData object containing b.
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', bt, 'h');
        end
        function set_normalisation(self, bin_eff)
%***SIRF*** Sets the normalisation n in (F);
%         bin_eff:  an AcquisitionData object containing bin efficiencies
%                   (the inverse of n).
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', bin_eff, 'h');
        end
        function set_up(self, acq_templ, img_templ)
%***SIRF*** Prepares this object for performing forward and backward
%         projections;
%         acq_templ:  an AcquisitionData object used as a template for
%                     creating an AcquisitionData object to store forward
%                     projection;
%         img_templ:  an ImageData object used as a template for creating an
%                     ImageData object to store backward projection.
            h = calllib...
                ('mstir', 'mSTIR_setupAcquisitionModel',...
                self.handle, acq_templ.handle, img_templ.handle);
            mUtil.checkExecutionStatus([self.name ':set_up'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function ad = forward(self, image, filename)
%***SIRF*** Returns the forward projection of x given by (F);
%         image   :  an ImageData object containing x;
%         filename:  an optional name of the file to store projection data;
%                    if not present, projection data is stored in memory
%                    (not recommended as it can be huge).
            if nargin < 3
                filename = '';
            end
            ad = mStir.AcquisitionData();
            ad.handle = calllib('mstir', 'mSTIR_acquisitionModelFwd',...
                self.handle, image.handle, filename);
            mUtil.checkExecutionStatus([self.name ':forward'], ad.handle)
        end
        function image = backward(self, ad)
%***SIRF*** Returns the backward projection of y giben by (B);
%         ad:  an AcquisitionData object containing y.
            image = mStir.ImageData();
            image.handle = calllib('mstir', 'mSTIR_acquisitionModelBwd',...
                self.handle, ad.handle);
            mUtil.checkExecutionStatus...
                ([self.name ':backward'], image.handle)
        end
    end
end