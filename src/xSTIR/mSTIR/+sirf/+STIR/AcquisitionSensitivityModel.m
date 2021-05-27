classdef AcquisitionSensitivityModel < handle
%     Class for PET acquisition sensitivity model objects.

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
        name_
        handle_
    end
    methods (Static)
        function name = class_name()
            name = 'AcquisitionSensitivityModel';
        end
    end
    methods
        function self = AcquisitionSensitivityModel(src, other_src)
%***SIRF*** Creates new AcquisitionSensitivityModel object
%         - from a manufacturer normalisation file (supported by STIR) or
%         - from ImageData object containing attenuation image (units: 1/cm) or
%         - from AcquisitionData object containing bin efficiencies or
%         - by chaining two existing AcquisitionSensitivityModel objects
%         src: file name or ImageData object or AcquisitionData object
%         other_src: AcquisitionSensitivityModel object (optional)
            self.handle_ = [];
            self.name_ = 'AcquisitionSensitivityModel';
            if nargin < 1
                return
            end
            if nargin < 2
                if ischar(src)
                    % create from ECAT8/GE norm file
                    fprintf(1, 'Reading manufacturer PET normalisation file from %s', src)
                    h = calllib('miutilities', 'mCharDataHandle', src);
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAcquisitionSensitivityModel', h, 'n');
                elseif isa(src, 'sirf.STIR.AcquisitionData')
                    assert(~isempty(src.handle_), 'empty bin efficiency data')
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAcquisitionSensitivityModel',...
                        src.handle_, 's');
                else
                    error([self.name_ ':wrong_ctor_source'], ...
                    'wrong source in AcquisitionSensitivityModel constructor')
                end
            else
                if isa(src, 'sirf.STIR.ImageData')
                    assert(~isempty(src.handle_), 'empty attenuation image')
                    assert(~isempty(other_src.handle_),...
                        'empty acquisition model')
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAttenuationModel',...
                        src.handle_, other_src.handle_);
                elseif isa(src, 'sirf.STIR.AcquisitionSensitivityModel') && ...
                        isa(other_src, 'sirf.STIR.AcquisitionSensitivityModel')
                    assert(~isempty(src.handle_),...
                        'first sensitivity model is empty')
                    assert(~isempty(other_src.handle_),...
                        'second sensitivity model is empty')
                    self.handle_ = calllib('mstir',...
                        'mSTIR_chainPETAcquisitionSensitivityModels',...
                        src.handle_, other_src.handle_);
                else
                    error([self.name_ ':wrong_ctor_source'], ...
                    'wrong sources in AcquisitionSensitivityModel constructor')
                end
            end
            sirf.Utilities.check_status([self.name_ ':ctor'], self.handle_)
        end
        function set_up(self, acq_data)
%***SIRF*** Sets up the object.
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
            h = calllib('mstir',...
                'mSTIR_setupAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_);
            sirf.Utilities.check_status([self.name_ ':set_up'], h)
            sirf.Utilities.delete(h)
        end
        function normalise(self, acq_data)
%***SIRF*** Multiplies the argument by n (cf. AcquisitionModel).
%         If self is a chain of two AcquisitionSensitivityModels, then 
%         n is a product of two normalisations.
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
            h = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'normalise');
            sirf.Utilities.check_status([self.name_ ':set_up'], h)
            sirf.Utilities.delete(h)
        end
        function unnormalise(self, acq_data)
%***SIRF*** Multiplies the argument by 1/n (cf. AcquisitionModel).
%         If self is a chain of two AcquisitionSensitivityModels, then 
%         n is a product of two normalisations.
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
            h = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'unnormalise');
            sirf.Utilities.check_status([self.name_ ':set_up'], h)
            sirf.Utilities.delete(h)
        end
        function fwd_data = forward(self, acq_data)
%***SIRF*** Same as unnormalise except that the argument remains unchanged
%         and a new AcquisitionData equal to the argument multiplied
%         by 1/n is returned.
           assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
            fwd_data = sirf.STIR.AcquisitionData();
            fwd_data.handle_ = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'fwd');
            sirf.Utilities.check_status([self.name_ ':set_up'], fwd_data.handle_)
        end
        function inv_data = invert(self, acq_data)
%***SIRF*** Same as normalise except that the argument remains unchanged
%         and a new AcquisitionData equal to the argument multiplied
%         by n is returned.
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            sirf.Utilities.assert_validity(acq_data, 'AcquisitionData')
            inv_data = sirf.STIR.AcquisitionData();
            inv_data.handle_ = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'inv');
            sirf.Utilities.check_status([self.name_ ':set_up'], inv_data.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
    end
end
