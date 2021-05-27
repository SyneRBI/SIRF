function check_status(f, handle)
% Checks the execution status recorded in handle.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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

    %lib = 'mutilities';
    lib = 'miutilities';
    status = calllib(lib, 'mExecutionStatus', handle);
    if status ~= 0
        errId = [f ':error'];
        msg = calllib(lib, 'mExecutionError', handle);
        file = calllib(lib, 'mExecutionErrorFile', handle);
        line = calllib(lib, 'mExecutionErrorLine', handle);
		msg2 = 'the reconstruction engine output may provide more information';
        error(errId, '??? %s exception thrown at line %d of %s, \n%s', ...
            msg, line, file, msg2)
    end
end
