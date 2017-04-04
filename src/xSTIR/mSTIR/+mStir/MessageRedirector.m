classdef MessageRedirector < handle
% Class for STIR error, warning and information messages redirection.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
        info
        info_case
        warn
        warn_case
        errr
        errr_case
    end
    methods
        function self = MessageRedirector(info, warn, errr)
%         MessageRedirector(info, warn, errr) creates a MessageRedirector object.
%         All arguments are either 'stdout' (default, messages are printed
%         in Matlab Command Window) or filenames.
%         info: destination for messages printed by STIR function info
%         warn: destination for messages printed by STIR function warning
%         errr: destination for messages printed by STIR function error
            if nargin > 0
                if isempty(info)
                    self.info_case = -1;
                elseif strcmpi(info, 'stdout') ~= 0
                    self.info = calllib('mstir', 'mNewMexPrinter');
                    self.info_case = 0;
                else
                    self.info = calllib('mstir', 'mNewTextWriter', info);
                    self.info_case = 1;
                end
                if self.info_case ~= -1
                    calllib('mstir', 'mOpenChannel', 0, self.info)
                else
                    calllib('mstir', 'mCloseChannel', 0, self.info)
                end
            else
                self.info = calllib('mstir', 'mNewMexPrinter');
                self.info_case = 0;
                calllib('mstir', 'mOpenChannel', 0, self.info)
            end
            if nargin > 1
                if isempty(warn)
                    self.warn_case = -1;
                    calllib('mstir', 'mCloseChannel', 1, self.warn)
                elseif strcmpi(warn, 'stdout') ~= 0
                    self.warn = calllib('mstir', 'mNewMexPrinter');
                    self.warn_case = 0;
                    calllib('mstir', 'mOpenChannel', 1, self.warn)
                else
                    self.warn = calllib('mstir', 'mNewTextWriter', warn);
                    self.warn_case = 1;
                    calllib('mstir', 'mOpenChannel', 1, self.warn)
                end
            else
                self.warn = calllib('mstir', 'mNewMexPrinter');
                self.warn_case = 0;
                calllib('mstir', 'mOpenChannel', 1, self.warn)
            end
            if nargin > 2
                if isempty(errr)
                    self.errr_case = -1;
                    calllib('mstir', 'mCloseChannel', 2, self.errr)
                elseif strcmpi(errr, 'stdout') ~= 0
                    self.errr = calllib('mstir', 'mNewMexPrinter');
                    self.errr_case = 0;
                    calllib('mstir', 'mOpenChannel', 2, self.errr)
                else
                    self.errr = calllib('mstir', 'mNewTextWriter', errr);
                    self.errr_case = 1;
                    calllib('mstir', 'mOpenChannel', 2, self.errr)
                end
            else
                self.errr = calllib('mstir', 'mNewMexPrinter');
                self.errr_case = 0;
                calllib('mstir', 'mOpenChannel', 2, self.errr)
            end
        end
        function delete(self)
            if self.info_case ~= -1
                calllib('mstir', 'mCloseChannel', 0, self.info)
                if self.info_case == 0
                    calllib('mstir', 'mDeleteMexPrinter', self.info)
                else
                    calllib('mstir', 'mDeleteTextWriter', self.info)
                end
            end
            if self.warn_case ~= -1
                calllib('mstir', 'mCloseChannel', 1, self.warn)
                if self.warn_case == 0
                    calllib('mstir', 'mDeleteMexPrinter', self.warn)
                else
                    calllib('mstir', 'mDeleteTextWriter', self.warn)
                end
            end
            if self.errr_case ~= -1
                calllib('mstir', 'mCloseChannel', 2, self.errr)
                if self.errr_case == 0
                    calllib('mstir', 'mDeleteMexPrinter', self.errr)
                else
                    calllib('mstir', 'mDeleteTextWriter', self.errr)
                end
            end
        end
    end
end