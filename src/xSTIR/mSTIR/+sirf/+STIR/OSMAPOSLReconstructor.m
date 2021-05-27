classdef OSMAPOSLReconstructor < sirf.STIR.IterativeReconstructor
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
        name
    end
    methods
        function self = OSMAPOSLReconstructor(filename)
%         Creates an OSMAPOSL reconstructor object.
%         The optional parameter provides the name of the Interfile setting
%         the parameters of reconstruction.
            self.name = 'OSMAPOSL';
            if nargin < 1
                filename = '';
            end
            self.handle_ = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'OSMAPOSLReconstruction', filename);
            sirf.Utilities.check_status(self.name, self.handle_);
        end
        function delete(self)
            %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
            sirf.Utilities.delete(self.handle_)
            self.handle_ = [];
        end
        function set_MAP_model(self, model)
            sirf.STIR.setParameter(self.handle_, self.name, 'MAP_model', model, 'c')
        end
    end
end
