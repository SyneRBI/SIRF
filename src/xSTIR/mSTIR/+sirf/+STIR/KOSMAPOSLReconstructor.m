classdef KOSMAPOSLReconstructor < sirf.STIR.IterativeReconstructor
%     Class for reconstructor objects using Kernel Ordered Subsets Maximum
%     A Posteriori One Step Late reconstruction algorithm.
% 
%     This class implements the iterative algorithm obtained using the Kernel
%     method (KEM) and Hybrid kernel method (HKEM). This implementation
%     corresponds to the one presented by Deidda D et al, "Hybrid PET-MR
%     list-mode kernelized expectation maximization  reconstruction", Inverse
%     Problems, 2019, DOI: https://doi.org/10.1088/1361-6420/ab013f.
%     However, this allows also sinogram-based reconstruction. Each voxel value
%     of the image X can be represented as a linear combination using the kernel
%     method.  If we have an image with prior information, we can construct for
%     each voxel j of the emission image a feature vector, v, using the prior
%     information. The image X can then be described using the kernel matrix
% 
%     X = A*K
% 
%     where K is the kernel matrix. The resulting algorithm with OSEM,
%     for example, is the following:
% 
%     A^(n+1) =  A^n/(K^n * S) * K^n * P * Y/(P * K^n *A^n + S)
% 
%     where kernel can be written as:
% 
%     K^n = K_m * K_p;
% 
%     with
% 
%     K_m = exp(-(v_j - v_l)^2/(2*sigma_m^2)) *
%           exp(-(x_j - x_l)^2 /(2*sigma_dm^2))
% 
%     being the MR component of the kernel and
% 
%     K_p = exp(-(z_j - z_l)^2/(2*sigma_p^2)) *
%           exp(-(x_j - x_l)^2 /(2*sigma_dp^2))
% 
%     is the part coming from the emission iterative update. Here, the Gaussian
%     kernel functions have been modulated by the distance between voxels in the
%     image space.

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
