classdef Image < handle
% INTERNAL USE ONLY.
% Class for the ISMRMRD image object.

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
            name = 'Image';
        end
        function obj = same_object()
            obj = mGadgetron.Image();
        end
    end
    methods
        function self = Image()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function v = version(self)
            assert(~isempty(self.handle_), 'empty Image object')
            v = mGadgetron.parameter(self.handle_, ...
                'image', 'version', 'i');
        end
        function f = flags(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'flags', 'i');
        end
        function f = data_type(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'data_type', 'i');
        end
        function uid = measurement_uid(self)
            assert(~isempty(self.handle_), 'empty Image object')
            uid = mGadgetron.parameter(self.handle_, ...
                'image', 'measurement_uid', 'i');
        end
        function nc = channels(self)
            assert(~isempty(self.handle_), 'empty Image object')
            nc = mGadgetron.parameter(self.handle_, ...
                'image', 'channels', 'i');
        end
        function a = average(self)
            assert(~isempty(self.handle_), 'empty Image object')
            a = mGadgetron.parameter(self.handle_, ...
                'image', 'average', 'i');
        end
        function s = slice(self)
            assert(~isempty(self.handle_), 'empty Image object')
            s = mGadgetron.parameter(self.handle_, ...
                'image', 'slice', 'i');
        end
        function c = contrast(self)
            assert(~isempty(self.handle_), 'empty Image object')
            c = mGadgetron.parameter(self.handle_, ...
                'image', 'contrast', 'i');
        end
        function p = phase(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'phase', 'i');
        end
        function r = repetition(self)
            assert(~isempty(self.handle_), 'empty Image object')
            r = mGadgetron.parameter(self.handle_, ...
                'image', 'repetition', 'i');
        end
        function s = set(self)
            assert(~isempty(self.handle_), 'empty Image object')
            s = mGadgetron.parameter(self.handle_, ...
                'image', 'set', 'i');
        end
        function f = image_type(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'image_type', 'i');
        end
        function f = image_index(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'image_index', 'i');
        end
        function f = image_series_index(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'image_series_index', 'i');
        end
        function f = attribute_string_len(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'attribute_string_len', 'i');
        end
        function f = matrix_size(self)
            assert(~isempty(self.handle_), 'empty Image object')
            f = mGadgetron.parameter(self.handle_, ...
                'image', 'matrix_size', 'u16s', 3);
        end
        function ats = acquisition_time_stamp(self)
            assert(~isempty(self.handle_), 'empty Image object')
            ats = mGadgetron.parameter(self.handle_, ...
                'image', 'acquisition_time_stamp', 'i');
        end
        function pts = physiology_time_stamp(self)
            assert(~isempty(self.handle_), 'empty Image object')
            pts = mGadgetron.parameter(self.handle_, ...
                'image', 'physiology_time_stamp', 'u32s', 3);
        end
        function p = field_of_view(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'field_of_view', 'fs', 3);
        end
        function p = position(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'position', 'fs', 3);
        end
        function p = read_dir(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'read_dir', 'fs', 3);
        end
        function p = phase_dir(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'phase_dir', 'fs', 3);
        end
        function p = slice_dir(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'slice_dir', 'fs', 3);
        end
        function p = patient_table_position(self)
            assert(~isempty(self.handle_), 'empty Image object')
            p = mGadgetron.parameter(self.handle_, ...
                'image', 'patient_table_position', 'fs', 3);
        end
        function i = info(self, method)
            i = eval(['self.' method '()']);
        end
    end
end