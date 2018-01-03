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
    methods (Static)
        function name = class_name()
            name = 'Acquisition';
        end
        function obj = same_object()
            obj = mGadgetron.Acquisition();
        end
    end
    methods
        function self = Acquisition()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function v = version(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            v = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'version', 'i');
        end
        function f = flags(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            f = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'flags', 'i');
        end
        function uid = measurement_uid(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            uid = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'measurement_uid', 'i');
        end
        function sc = scan_counter(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            sc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'scan_counter', 'i');
        end
        function ats = acquisition_time_stamp(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            ats = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'acquisition_time_stamp', 'i');
        end
        function ns = number_of_samples(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'number_of_samples', 'i');
        end
        function ac = available_channels(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            ac = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'available_channels', 'i');
        end
        function nc = active_channels(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            nc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'active_channels', 'i');
        end
        function dp = discard_pre(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            dp = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'discard_pre', 'i');
        end
        function dp = discard_post(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            dp = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'discard_post', 'i');
        end
        function cs = center_sample(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            cs = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'center_sample', 'i');
        end
        function ref = encoding_space_ref(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            ref = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'encoding_space_ref', 'i');
        end
        function td = trajectory_dimensions(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            td = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'trajectory_dimensions', 'i');
        end
        function es = kspace_encode_step_1(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_1', 'i');
        end
        function es = kspace_encode_step_2(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_2', 'i');
        end
        function a = average(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            a = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_average', 'i');
        end
        function s = slice(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_slice', 'i');
        end
        function c = contrast(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            c = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_contrast', 'i');
        end
        function p = phase(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_phase', 'i');
        end
        function r = repetition(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_repetition', 'i');
        end
        function s = set(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_set', 'i');
        end
        function s = segment(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_segment', 'i');
        end
        function pts = physiology_time_stamp(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            pts = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'physiology_time_stamp', 'u32s', 3);
        end
        function m = channel_mask(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            m = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'channel_mask', 'u64s', 16);
        end
        function s = sample_time_us(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            s = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'sample_time_us', 'f');
        end
        function p = position(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'position', 'fs', 3);
        end
        function p = read_dir(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'read_dir', 'fs', 3);
        end
        function p = phase_dir(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'phase_dir', 'fs', 3);
        end
        function p = slice_dir(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'slice_dir', 'fs', 3);
        end
        function p = patient_table_position(self)
            assert(~isempty(self.handle_), 'empty Acquisition object')
            p = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'patient_table_position', 'fs', 3);
        end
        function i = info(self, method)
            i = eval(['self.' method '()']);
        end
    end
end