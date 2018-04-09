classdef PLSPrior < mSTIR.Prior
% Class for PLS priors

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

    methods
        function self = PLSPrior()
            self.name = 'PLSPrior';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
            mUtilities.check_status('PLSPrior:ctor', self.handle_)
        end
        function delete(self)
            mUtilities.delete(self.handle_)
            self.handle_ = [];
        end
        function set_only_2D(self, name)
%***SIRF*** Sets only 2D.
            mSTIR.setParameter(self.handle_, 'PLSPrior', 'only_2D', name, 'f')
        end
        function set_alpha(self, name)
%***SIRF*** Sets alpha.
            mSTIR.setParameter(self.handle_, 'PLSPrior', 'alpha', name, 'f')
        end
        function set_eta(self, name)
%***SIRF*** Sets eta.
            mSTIR.setParameter(self.handle_, 'PLSPrior', 'eta', name, 'f')
        end
        function set_kappa_filename(self, name)
%***SIRF*** Sets kappa filename.
            mSTIR.setParameter(self.handle_, 'PLSPrior', 'kappa_filename', name, 'c')
        end
        function set_anatomical_filename(self, name)
%***SIRF*** Sets anatomical filename.
            mSTIR.setParameter(self.handle_, 'PLSPrior', 'anatomical_filename', name, 'c')
        end
        function value = get_only_2D(self)
%***SIRF*** Gets only 2D.
            value = mSTIR.parameter(self.handle_, 'PLSPrior', 'only_2D', 'f');
        end
        function value = get_alpha(self)
%***SIRF*** Gets alpha.
            value = mSTIR.parameter(self.handle_, 'PLSPrior', 'alpha', 'f');
        end
        function value = get_eta(self)
%***SIRF*** Gets eta.
            value = mSTIR.parameter(self.handle_, 'PLSPrior', 'eta', 'f');
        end
        function value = get_kappa_filename(self)
%***SIRF*** Gets kappa filename.
            value = mSTIR.parameter(self.handle_, 'PLSPrior', 'kappa_filename', 'c');
        end
        function value = get_anatomical_filename(self)
%***SIRF*** Gets anatomical filename.
            value = mSTIR.parameter(self.handle_, 'PLSPrior', 'anatomical_filename', 'c');
        end
        %function set_up(self)
%***SIRF*** Set up the PLS prior.
        %    h = calllib('mstir', 'mSTIR_setupPLSPrior', ...
        %        self.handle_);
        %    mUtilities.check_status([self.name_ ':set_up'], h);
        %    mUtilities.delete(h)
        %end
    end
end