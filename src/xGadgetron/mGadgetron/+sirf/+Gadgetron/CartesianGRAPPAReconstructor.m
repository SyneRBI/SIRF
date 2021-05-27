classdef CartesianGRAPPAReconstructor < sirf.Gadgetron.Reconstructor
% Class for the GRAPPA reconstructor from undersampled Cartesian raw data.

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
    end
    methods
        function self = CartesianGRAPPAReconstructor()
%         Creates Cartesian GRAPPA Reconstructor.
            self.name_ = 'SimpleGRAPPAReconstructionProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
            sirf.Utilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function compute_gfactors(self, flag)
%***SIRF*** Switches the computation of G-factors during the reconstruction
%         on (flag = true) or off (flag = false).
            self.set_gadget_property('gadget4', 'send_out_gfactor', flag)
        end
    end
end