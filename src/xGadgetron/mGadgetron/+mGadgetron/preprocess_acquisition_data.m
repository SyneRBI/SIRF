function prep_data = preprocess_acquisition_data(input_data)
% Pre-process the acquisition data specified by the argument.
% The following pre-processing steps are carried out:
%   - Noise prewhitening (NoiseAdjustGadget)
%   - Compensate for asymmetric echo acquisition (AsymmetricEchoAdjustROGadget) 
%   - Remove oversampling along acquisition direction (RemoveROOversamplingGadget)


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

prep_data = input_data.process([ ...
    {'NoiseAdjustGadget'} ...
    {'AsymmetricEchoAdjustROGadget'} ...
    {'RemoveROOversamplingGadget'} ...
    ]);
end
