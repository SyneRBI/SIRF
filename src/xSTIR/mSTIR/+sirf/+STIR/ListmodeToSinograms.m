classdef ListmodeToSinograms < handle
% Class for a listmode-to-sinograms converter
% This class reads list mode data and produces corresponding *sinograms*,
% i.e. histogrammed data in the format of PETAcquisitionData.
% It has two main functions:
%   - process() can be used to read prompts and/or delayed coincidences to
%     produce a single PETAcquisitionData.
%     Two conversion flags decide what is to be done with 3 possible cases:
%     - `store_prompts`=`true`, `store_delayeds`=`false`: only prompts stored
%     - `store_prompts`=`false`, `store_delayeds`=`true`: only delayeds stored
%     - `store_prompts`=`true`, `store_delayeds`=`true`: prompts-delayeds stored
%     Clearly, enabling the `store_delayeds` option only makes sense if the
%     data was acquired accordingly.
%   - estimate_randoms() can be used to get a relatively noiseless estimate of the 
%     random coincidences. 
% Currently, the randoms are estimated from the delayed coincidences using the
% following strategy:
%    1. singles (one per detector) are estimated using a Maximum Likelihood
%       estimator
%    2. randoms-from-singles are computed per detector-pair via the usual
%       product formula. These are then added together for all detector pairs
%       in a certain histogram-bin in the data (accommodating for view mashing
%       and axial compression).
% 
% The actual algorithm is described in
% D. Hogg, K. Thielemans, S. Mustafovic, and T. J. Spinks,
% "A study of bias for various iterative reconstruction methods in PET,"
% in 2002 IEEE Nuclear Science Symposium Conference Record, vol. 3. IEEE,
% Nov. 2002, pp. 1519-1523 (http://dx.doi.org/10.1109/nssmic.2002.1239610).

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
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_input(self, file)
            %***SIRF*** Sets the listmode file name.
            sirf.STIR.setParameter(self.handle_, self.name_, 'input', file, 'c')
        end
        function set_output_prefix(self, file)
            %***SIRF*** Sets the sinograms file names prefix.
            sirf.STIR.setParameter(self.handle_, self.name_, 'output', file, 'c')
        end
        function set_template(self, file)
            %***SIRF*** Sets the sinograms template.
            sirf.STIR.setParameter(self.handle_, self.name_, 'template', file, 'c')
        end
        function set_time_interval(self, start, stop)
            %***SIRF*** Sets time interval.
            % Only data scanned during this time interval will be converted.
            ptr = libpointer('singlePtr', [start stop]);
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsInterval', ...
                self.handle_, ptr);
            sirf.Utilities.check_status([self.name_ ':set_interval'], h);
            sirf.Utilities.delete(h)
        end
        function flag_on(self, flag)
            %***SIRF*** Switches on (sets to 'true') a conversion flag 
            % (see conversion flags description above).
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsFlag', ...
                self.handle_, flag, 1);
            sirf.Utilities.check_status([self.name_ ':flag_on'], h);
            sirf.Utilities.delete(h)
        end
        function flag_off(self, flag)
            %***SIRF*** Switches off (sets to 'false') a conversion flag 
            % (see conversion flags description above).
            h = calllib('mstir', 'mSTIR_setListmodeToSinogramsFlag', ...
                self.handle_, flag, 0);
            sirf.Utilities.check_status([self.name_ ':flag_on'], h);
            sirf.Utilities.delete(h)
        end
        function set_up(self)
            %***SIRF*** Sets up the conversion.
            h = calllib('mstir', 'mSTIR_setupListmodeToSinogramsConverter', ...
                self.handle_);
            sirf.Utilities.check_status([self.name_ ':set_up'], h);
            sirf.Utilities.delete(h)
        end
        function process(self)
            %***SIRF*** Performs the conversion.
            self.output_ = sirf.STIR.AcquisitionData();
            self.output_.handle_ = calllib...
                ('mstir', 'mSTIR_convertListmodeToSinograms', ...
                self.handle_);
            sirf.Utilities.check_status...
                ([self.name_ ':process'], self.output_.handle_);
        end
        function output = get_output(self)
            %***SIRF*** Returns the sinograms.
            assert(~isempty(self.output_), 'Conversion to sinograms not done')
            output = self.output_;
        end
        function randoms = estimate_randoms(self)
            %***SIRF*** Estimates randoms.
            randoms = sirf.STIR.AcquisitionData();
            randoms.handle_ = calllib('mstir', 'mSTIR_computeRandoms', ...
                self.handle_);
            sirf.Utilities.check_status...
                ([self.name_ ':estimate_randoms'], randoms.handle_);
        end
        function v = get_time_at_which_prompt_rate_exceeds_threshold(self, threshold)
            %Get the time in the list mode data at which the number
            %of prompts per second exceeds a given threshold.
            %Returns -1 if no corresponding time is found.
            h = calllib('mstir', 'mSTIR_lm_prompt_rate_exceeds_threshold',...
                            self.handle_, threshold);
            sirf.Utilities.check_status...
                ([self.name_ '::get_time_at_which_prompt_rate_exceeds_threshold'], h)
            v = calllib('miutilities', 'mFloatDataFromHandle', h);
            sirf.Utilities.delete(h)
        end
    end
end