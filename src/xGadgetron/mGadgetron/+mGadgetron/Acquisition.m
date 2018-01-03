classdef Acquisition < handle
% INTERNAL USE ONLY.
% Class for the ISMRMRD acquisition object.

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
    end
    methods
        function self = Acquisition()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function v = version(self)
            v = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'version', 'i');
        end
        function f = flags(self)
            f = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'flags', 'i');
        end
        function uid = measurement_uid(self)
            uid = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'measurement_uid', 'i');
        end
        function sc = scan_counter(self)
            sc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'scan_counter', 'i');
        end
        function ats = acquisition_time_stamp(self)
            ats = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'acquisition_time_stamp', 'i');
        end
        function ns = number_of_samples(self)
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'number_of_samples', 'i');
        end
        function ac = available_channels(self)
            ac = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'available_channels', 'i');
        end
        function nc = active_channels(self)
            nc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'active_channels', 'i');
        end
        function dp = discard_pre(self)
            dp = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'discard_pre', 'i');
        end
        function dp = discard_post(self)
            dp = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'discard_post', 'i');
        end
        function cs = center_sample(self)
            cs = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'center_sample', 'i');
        end
        function ref = encoding_space_ref(self)
            ref = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'encoding_space_ref', 'i');
        end
        function td = trajectory_dimensions(self)
            td = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'trajectory_dimensions', 'i');
        end
        function es = kspace_encode_step_1(self)
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_1', 'i');
        end
        function es = kspace_encode_step_2(self)
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_2', 'i');
        end
        function a = average(self)
            a = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_average', 'i');
        end
        function s = slice(self)
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_slice', 'i');
        end
        function c = contrast(self)
            c = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_contrast', 'i');
        end
        function p = phase(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_phase', 'i');
        end
        function r = repetition(self)
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_repetition', 'i');
        end
        function s = set(self)
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_set', 'i');
        end
        function s = segment(self)
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_segment', 'i');
        end
        function pts = physiology_time_stamp(self)
            pts = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'physiology_time_stamp', 'u32s', 3);
        end
        function m = channel_mask(self)
            m = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'channel_mask', 'u64s', 16);
        end
        function s = sample_time_us(self)
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'sample_time_us', 'f');
        end
        function p = position(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'position', 'fs', 3);
        end
        function p = read_dir(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'read_dir', 'fs', 3);
        end
        function p = phase_dir(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'phase_dir', 'fs', 3);
        end
        function p = slice_dir(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'slice_dir', 'fs', 3);
        end
        function p = patient_table_position(self)
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'patient_table_position', 'fs', 3);
        end
        function i = info(self, method)
            i = eval(['self.' method]);
        end
    end
end