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
        output_
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
            % Sets the listmode file name.
            mSTIR.setParameter(self.handle_, self.name_, 'input', file, 'c')
        end
        function set_output_prefix(self, file)
            % Sets the sinograms file names prefix.
            mSTIR.setParameter(self.handle_, self.name_, 'output', file, 'c')
        end
        function set_template(self, file)
            % Sets the sinograms template.
            mSTIR.setParameter(self.handle_, self.name_, 'template', file, 'c')
        end
        function set_time_interval(self, start, stop)
            % Sets time interval.
            % Only data scanned during this time interval will be converted.
            ptr = libpointer('singlePtr', [start stop]);
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsInterval', ...
                self.handle_, ptr);
            mUtilities.check_status([self.name_ ':set_interval'], h);
            mUtilities.delete(h)
        end
        function flag_on(self, flag)
            % Switches on a conversion flag.
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsFlag', ...
                self.handle_, flag, 1);
            mUtilities.check_status([self.name_ ':flag_on'], h);
            mUtilities.delete(h)
        end
        function flag_off(self, flag)
            % Switches off a conversion flag.
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsFlag', ...
                self.handle_, flag, 0);
            mUtilities.check_status([self.name_ ':flag_on'], h);
            mUtilities.delete(h)
        end
        function set_up(self)
            % Sets up the conversion.
            h = calllib('mstir', 'mSTIR_setupListmodeToSinogramsConverter', ...
                self.handle_);
            mUtilities.check_status([self.name_ ':set_up'], h);
            mUtilities.delete(h)
        end
        function process(self)
            % Performs the conversion.
            self.output_ = mSTIR.AcquisitionData();
            self.output_.handle_ = calllib...
                ('mstir', 'mSTIR_convertListmodeToSinograms', ...
                self.handle_);
            mUtilities.check_status...
                ([self.name_ ':process'], self.output_.handle_);
        end
        function output = get_output(self)
            % Returns the sinograms.
            assert(~isempty(self.output_), 'Conversion to sinograms not done')
            output = self.output_;
        end
        function randoms = estimate_randoms(self)
            % Estimates randoms.
            randoms = mSTIR.AcquisitionData();
            randoms.handle_ = calllib('mstir', 'mSTIR_computeRandoms', ...
                self.handle_);
            mUtilities.check_status...
                ([self.name_ ':estimate_randoms'], randoms.handle_);
        end
    end
end