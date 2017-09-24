function assert_validities(x, y)
% Ensures x and y are of the same type and not empty.

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

assert(~isempty(x.handle_), 'first object is empty')
assert(~isempty(y.handle_), 'second object is empty')
if ~strcmp(class(x), class(y))
    fprintf('??? Objects types are %s and %s: same type expected.\n', ...
        class(x), class(y))
    error('Objects must be of the same type')
end
%assert(strcmp(class(x), class(y)), 'objects must be of the same type')