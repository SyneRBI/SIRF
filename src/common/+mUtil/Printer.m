classdef Printer < handle
    properties
        info
        info_case
        warn
        warn_case
        errr
        errr_case
    end
    methods
        function self = Printer(info, warn, errr)
            lib = 'mutilities';
            if nargin > 0
                if isempty(info)
                    self.info_case = -1;
                elseif strcmpi(info, 'stdout') ~= 0
                    self.info = calllib(lib, 'mNewMexPrinter');
                    self.info_case = 0;
                else
                    self.info = calllib(lib, 'mNewTextWriter', info);
                    self.info_case = 1;
                end
                if self.info_case ~= -1
                    calllib(lib, 'mOpenChannel', 0, self.info)
                else
                    calllib(lib, 'mCloseChannel', 0, self.info)
                end
            end
            if nargin > 1
                if isempty(warn)
                    self.warn_case = -1;
                    calllib(lib, 'mCloseChannel', 1, self.warn)
                elseif strcmpi(warn, 'stdout') ~= 0
                    self.warn = calllib(lib, 'mNewMexPrinter');
                    self.warn_case = 0;
                    calllib(lib, 'mOpenChannel', 1, self.warn)
                else
                    self.warn = calllib(lib, 'mNewTextWriter', warn);
                    self.warn_case = 1;
                    calllib(lib, 'mOpenChannel', 1, self.warn)
                end
            else
                self.warn = calllib(lib, 'mNewMexPrinter');
                self.warn_case = 0;
                calllib(lib, 'mOpenChannel', 1, self.warn)
            end
            if nargin > 2
                if isempty(errr)
                    self.errr_case = -1;
                    calllib(lib, 'mCloseChannel', 2, self.errr)
                elseif strcmpi(errr, 'stdout') ~= 0
                    self.errr = calllib(lib, 'mNewMexPrinter');
                    self.errr_case = 0;
                    calllib(lib, 'mOpenChannel', 2, self.errr)
                else
                    self.errr = calllib(lib, 'mNewTextWriter', errr);
                    self.errr_case = 1;
                    calllib(lib, 'mOpenChannel', 2, self.errr)
                end
            else
                self.errr = calllib(lib, 'mNewMexPrinter');
                self.errr_case = 0;
                calllib(lib, 'mOpenChannel', 2, self.errr)
            end
        end
        function delete(self)
            lib = 'mutilities';
            if self.info_case ~= -1
                calllib(lib, 'mCloseChannel', 0, self.info)
                if self.info_case == 0
                    calllib(lib, 'mDeleteMexPrinter', self.info)
                else
                    calllib(lib, 'mDeleteTextWriter', self.info)
                end
            end
            if self.warn_case ~= -1
                calllib(lib, 'mCloseChannel', 1, self.warn)
                if self.warn_case == 0
                    calllib(lib, 'mDeleteMexPrinter', self.warn)
                else
                    calllib(lib, 'mDeleteTextWriter', self.warn)
                end
            end
            if self.errr_case ~= -1
                calllib(lib, 'mCloseChannel', 2, self.errr)
                if self.errr_case == 0
                    calllib(lib, 'mDeleteMexPrinter', self.errr)
                else
                    calllib(lib, 'mDeleteTextWriter', self.errr)
                end
            end
        end
    end
end