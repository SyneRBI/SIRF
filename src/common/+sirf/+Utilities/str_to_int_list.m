function int_list = str_to_int_list(str_list)
% Converts the string str_list of the form n1{-n2] [, n2[-n3]] [, ...]
% into a list of numbers n1[, n1 + 1, ..., n2] etc.

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

int_list = [];
last = false;
while ~last
    ic = strfind(str_list, ',');
    if isempty(ic)
        ic = length(str_list) + 1;
        last = true;
    end
    str_item = str_list(1 : ic - 1);
    str_list = str_list(ic + 1 : end);
    ic = strfind(str_item, '-');
    if isempty(ic)
        int_item = [str2num(str_item)];
    else
        strt = [str2num(str_item(1 : ic - 1))];
        stop = [str2num(str_item(ic + 1 : end))];
        int_item = strt : stop;
    end
    int_list = [int_list int_item];    
end
end