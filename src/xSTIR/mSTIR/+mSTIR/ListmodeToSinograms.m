classdef ListmodeToSinograms < handle
% Class for a listmode-to-sinograms converter

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
        handle_
        name_
    end
    methods
        function self = ListmodeToSinograms(file)
            self.name_ = 'ListmodeToSinograms';
            if nargin < 1
                self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name_);
            else
                self.handle_ = ...
                    calllib('mstir', 'mSTIR_newObject', self.name_, file);
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_input(self, file)
            mSTIR.setParameter(self.handle_, self.name_, 'input', file, 'c')
        end
        function set_output(self, file)
            mSTIR.setParameter(self.handle_, self.name_, 'output', file, 'c')
        end
        function set_template(self, file)
            mSTIR.setParameter(self.handle_, self.name_, 'template', file, 'c')
        end
        function set_interval(self, start, stop)
            ptr = libpointer('singlePtr', [start stop]);
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsInterval', ...
                self.handle_, ptr);
            mUtilities.check_status([self.name_ ':set_interval'], h);
            mUtilities.delete(h)
        end
        function set_up(self)
            h = calllib('mstir', 'mSTIR_setupListmodeToSinogramsConverter', ...
                self.handle_);
            mUtilities.check_status([self.name_ ':set_up'], h);
            mUtilities.delete(h)
        end
        function process(self)
            h = calllib('mstir', 'mSTIR_convertListmodeToSinograms', ...
                self.handle_);
            mUtilities.check_status([self.name_ ':set_up'], h);
            mUtilities.delete(h)
        end
    end
end