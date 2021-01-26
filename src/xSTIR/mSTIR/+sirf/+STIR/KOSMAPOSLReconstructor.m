classdef KOSMAPOSLReconstructor < sirf.STIR.IterativeReconstructor
% Class for reconstructor objects using OSMAPOSL algorithm.
% OSMAPOSL stands for Ordered Subsets Maximum A Posteriori One Step Late, see
% http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html

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
    end
    methods
        function self = KOSMAPOSLReconstructor(filename)
%         Creates an KOSMAPOSL reconstructor object.
%         The optional parameter provides the name of the Interfile setting
%         the parameters of reconstruction.
            self.name_ = 'KOSMAPOSL';
            if nargin < 1
                filename = '';
            end
            self.handle_ = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'KOSMAPOSLReconstruction', filename);
            sirf.Utilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            sirf.Utilities.delete(self.handle_)
            self.handle_ = [];
        end
        function set_anatomical_prior(self, ap)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'anatomical_prior', ap, 'h');
        end
        function set_num_neighbours(self, n)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'num_neighbours', n, 'i');
        end
        function set_num_non_zero_features(self, n)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'num_non_zero_features', n, 'i');
        end
        function set_sigma_m(self, v)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'sigma_m', v, 'f');
        end
        function set_sigma_p(self, v)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'sigma_p', v, 'f');
        end
        function set_sigma_dm(self, v)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'sigma_dm', v, 'f');
        end
        function set_sigma_dp(self, v)
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'sigma_dp', v, 'f');
        end
        function set_only_2D(self, tf)
            if tf
                v = 1;
            else
                v = 0;
            end
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'only_2D', v, 'i');
        end
        function set_hybrid(self, tf)
            if tf
                v = 1;
            else
                v = 0;
            end
            sirf.STIR.setParameter...
                (self.handle_, self.name_, 'hybrid', v, 'i');
        end
    end
end
