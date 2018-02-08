classdef AcquisitionSensitivityModel < handle
%     Class for PET acquisition sensitivity model objects.

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
%         Creates new AcquisitionSensitivityModel object
%         - from an ECAT8 file or
%         - from ImageData object containing attenuation image or
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
                    h = calllib('miutilities', 'mCharDataHandle', src);
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAcquisitionSensitivityModel', h, 'n');
                elseif isa(src, 'mSTIR.AcquisitionData')
                    assert(~isempty(src.handle_), 'empty bin efficiency data')
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAcquisitionSensitivityModel',...
                        src.handle_, 's');
                else
                    error([self.name_ ':wrong_ctor_source'], ...
                    'wrong source in AcquisitionSensitivityModel constructor')
                end
            else
                if isa(src, 'mSTIR.ImageData')
                    assert(~isempty(src.handle_), 'empty attenuation image')
                    assert(~isempty(other_src.handle_),...
                        'empty acquisition model')
                    self.handle_ = calllib('mstir',...
                        'mSTIR_createPETAttenuationModel',...
                        src.handle_, other_src.handle_);
                elseif isa(src, 'mSTIR.AcquisitionSensitivityModel') && ...
                        isa(other_src, 'mSTIR.AcquisitionSensitivityModel')
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
            mUtilities.check_status([self.name_ ':ctor'], self.handle_)
        end
        function set_up(self, acq_data)
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            mUtilities.assert_validity(acq_data, 'AcquisitionData')
            h = calllib('mstir',...
                'mSTIR_setupAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_);
            mUtilities.check_status([self.name_ ':set_up'], h)
            mUtilities.delete(h)
        end
        function unnormalise(self, acq_data)
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            mUtilities.assert_validity(acq_data, 'AcquisitionData')
            h = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'unnormalise');
            mUtilities.check_status([self.name_ ':set_up'], h)
            mUtilities.delete(h)
        end
        function fwd_data = forward(self, acq_data)
            assert(~isempty(self.handle_),...
                'empty acquisition sensitivity object')
            mUtilities.assert_validity(acq_data, 'AcquisitionData')
            fwd_data = mSTIR.AcquisitionData();
            fwd_data.handle_ = calllib('mstir',...
                'mSTIR_applyAcquisitionSensitivityModel',...
                self.handle_, acq_data.handle_, 'fwd');
            mUtilities.check_status([self.name_ ':set_up'], fwd_data.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
    end
end