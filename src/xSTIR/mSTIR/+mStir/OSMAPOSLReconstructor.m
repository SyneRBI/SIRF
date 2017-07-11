classdef OSMAPOSLReconstructor < mStir.IterativeReconstructor
% Class for reconstructor objects using OSMAPOSL algorithm.
% OSMAPOSL stands for Ordered Subsets Maximum A Posteriori One Step Late, see
% http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html

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
            self.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'OSMAPOSLReconstruction', filename);
            mUtil.check_status(self.name, self.handle);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
            self.handle = [];
        end
        function set_MAP_model(self, model)
            mStir.setParameter(self.handle, self.name, 'MAP_model', model, 'c')
        end
    end
end
