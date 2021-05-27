classdef CoilSensitivityData < sirf.SIRF.DataContainer
% Class for a coil sensitivity maps (csm) container.
% Each item in the container is a 4D (x-y-z-coil) complex array of csm 
% values.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
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
    end
    methods (Static)
        function name = class_name()
            name = 'CoilSensitivityData';
        end
        function obj = same_object()
            obj = sirf.Gadgetron.CoilSensitivityData();
        end
    end
    methods
        function self = CoilSensitivityData()
%         Creates an empty object.
            self.name_ = 'CoilSensitivityData';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function calculate(self, acqs)
%***SIRF*** Calculates coil sensitivity maps from sorted acquisitions
%         specified by an AcquisitionData argument.
            sirf.Utilities.assert_validity(acqs, 'AcquisitionData')
            if ~acqs.is_sorted()
                fprintf('WARNING: acquisitions may be in a wrong order\n')
            end
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = calllib('mgadgetron', 'mGT_CoilSensitivities', '');
            sirf.Utilities.check_status(self.name_, self.handle_);
            handle = calllib('mgadgetron', 'mGT_computeCoilSensitivities', ...
                self.handle_, acqs.handle_);
            sirf.Utilities.check_status(self.name_, handle);
            sirf.Utilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function data = csm_as_array(self, csm_num)
%***SIRF*** Returns the coil sensitivity map specified by the argument 
%         as a 4D (x-y-z-coil) complex array.
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib...
                ('mgadgetron', 'mGT_getCoilDataDimensions', ...
                self.handle_, csm_num - 1, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3)*dim(4);
            ptr_re = libpointer('singlePtr', zeros(n, 1));
            ptr_im = libpointer('singlePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getCoilData', ...
                self.handle_, csm_num - 1, ptr_re, ptr_im)
            re = reshape(ptr_re.Value, dim(1), dim(2), dim(3), dim(4));
            im = reshape(ptr_im.Value, dim(1), dim(2), dim(3), dim(4));
            data = complex(re, im);
        end
    end
end
