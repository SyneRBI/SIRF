classdef FBP2DReconstructor < sirf.STIR.Reconstructor
% Class for reconstructor objects using FBP2D algorithm.

% FBP2D stands for Filtered Back Projection 2D, see
% http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1FBP2DReconstruction.html

% Oblique angles in data will be ignored. The exception is the span=1 case,
% where the ring differences +1 and -1 are first combined to give indirect
% sinograms.
% By default, the algorithm uses the ramp filter. An apodizing filter can be
% added by using set_alpha_cosine_window and/or set_frequency_cut_off.
% The apodizing filter in frequency space has the form
% 
%     (alpha + (1 - alpha) * cos(pi * f / fc))

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
        function self = FBP2DReconstructor()
%         Creates an FBP2D reconstructor object.
%         The optional parameter provides the name of the Interfile setting
%         the parameters of reconstruction.
            self.name = 'FBP2D';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', 'FBP2D');
            sirf.Utilities.check_status(self.name, self.handle_);
        end
        function delete(self)
            sirf.Utilities.delete(self.handle_)
            self.handle_ = [];
        end
        function set_input(self, input_data)
            sirf.STIR.setParameter(self.handle_, self.name, 'input', input_data, 'h')
        end
        function set_zoom(self, zoom)
            sirf.STIR.setParameter(self.handle_, self.name, 'zoom', zoom, 'f')
        end
        function set_alpha_cosine_window(self, alpha)
%***SIRF*** Set alpha in the apodizing filter.
%         See the class documentation for the filter. The value of alpha should
%         be between 0.5 and 1. alpha=0.5 corresponds to the Hann filter, while
%         0.54 corresponds to the Hamming filter.
            sirf.STIR.setParameter(self.handle_, self.name, 'alpha', alpha, 'f')
        end
        function set_frequency_cut_off(self, fc)
%***SIRF*** Set the cut-off frequency for the apodizing filter.
%         See the class documentation for the filter. The value of fc should be
%         between 0 and 0.5.
            sirf.STIR.setParameter(self.handle_, self.name, 'fc', fc, 'f')
        end
        function set_output_image_size_xy(self, xy)
            sirf.STIR.setParameter(self.handle_, self.name, 'xy', xy, 'i')
        end
        function set_up(self, image)
            h = calllib('mstir', 'mSTIR_setupFBP2DReconstruction', ...
                self.handle_, image.handle_);
            sirf.Utilities.check_status([self.name ':set_up'], h);
        end
        function reconstruct(self)
            h = calllib('mstir', 'mSTIR_runFBP2DReconstruction', ...
                self.handle_);
            sirf.Utilities.check_status([self.name ':reconstruct'], h);
        end
        function image = get_output(self)
            image = sirf.STIR.ImageData();
            image.handle_ = calllib('mstir', 'mSTIR_parameter', ...
                self.handle_, self.name, 'output');
            sirf.Utilities.check_status([self.name ':get_output'], image.handle_);
        end
    end
end
