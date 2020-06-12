classdef mTest < handle
% Convenience class for writing test scripts.
    
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

    properties
        fid_
        data_
        record_
        rel_tol
        ntest
        failed
        verbose
    end
    methods
        function self = mTest(filename, record)
            self.fid_ = [];
            self.data_ = [];
            self.record_ = record;
            self.rel_tol = 1e-4;
            self.ntest = 0;
            self.failed = 0;
            self.verbose = true;
            if record
                self.fid_ = fopen(filename, 'w');
            else
                try
                    fid = fopen(filename, 'r');
                    self.data_ = fscanf(fid, '%e');
                    fclose(fid);
                catch
                    error('Error reading file %s\n', filename)
                end
            end
        end
        function delete(self)
            if ~isempty(self.fid_)
                fclose(self.fid_);
            end
        end
        function check(self, value, abs_tol)
            if nargin < 3
                abs_tol = 0.0;
            end
            self.ntest = self.ntest + 1;
            if self.record_
                fprintf(self.fid_, '%e\n', value);
            else
                if self.ntest > size(self.data_, 1)
                    fprintf('no data available for test %d\n', self.ntest);
                else
                    expected = self.data_(self.ntest);
                    eps = abs_tol + self.rel_tol*abs(expected);
                    if abs(value - expected) <= eps
                        if self.verbose
                            fprintf('+++ test %d passed\n', self.ntest);
                        end
                    else
                        self.failed = self.failed + 1;
                        if self.verbose
                            fprintf...
                                ('+++ test %d failed: expected %e, got %e\n', ...
                                self.ntest, expected, value);
                        end
                    end
                end
            end
        end
        function check_if_equal(self, expected, value)
            if value ~= expected
                self.failed = self.failed + 1;
                if self.verbose
                    fprintf...
                        ('+++ test %d failed: expected %e, got %e\n', ...
                        self.ntest, expected, value);
                end
            else
                if self.verbose
                    fprintf('+++ test %d passed\n', self.ntest);
                end
            end
            self.ntest = self.ntest + 1;
        end
    end
end