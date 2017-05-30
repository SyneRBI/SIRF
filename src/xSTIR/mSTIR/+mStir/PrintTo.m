classdef PrintTo < handle
    properties
        info
        info_case
        warn
        warn_case
        errr
        errr_case
    end
    methods
        function self = PrintTo(info, warn, errr)
%         PrintTo(info, warn, errr) redirects printing destinations.
%         info: specifies destination for information messages;
%         warn: specifies destination for warning messages;
%         errr: specifies destination for error messages.
%         All arguments can be:
%           'stdout' : printing goes to Matlab Command Window
%           string   : printing goes to a new file named by the string
%           empty    : printing is suppressed
            if nargin < 3 || isempty(errr)
                errr = '';
            end
            if nargin < 2 || isempty(warn)
                warn = '';
            end
            if nargin < 1 || isempty(info)
                info = '';
            end
            if nargin > 0
                if ~ischar(info)
                    error('Printer:wrong_arg', '??? Wrong info in Printer');
                elseif strcmpi(info, 'stdout') ~= 0
                    self.info = calllib('mstir', 'mNewMexPrinter');
                    self.info_case = 0;
                else
                    self.info = calllib('mstir', 'mNewTextWriter', info);
                    self.info_case = 1;
                end
                calllib('mstir', 'mOpenChannel', 0, self.info)
            end
            if nargin > 1
                if ~ischar(warn)
                    error('Printer:wrong_arg', '??? Wrong warn in Printer');
                elseif strcmpi(warn, 'stdout') ~= 0
                    self.warn = calllib('mstir', 'mNewMexPrinter');
                    self.warn_case = 0;
                else
                    self.warn = calllib('mstir', 'mNewTextWriter', warn);
                    self.warn_case = 1;
                end
                calllib('mstir', 'mOpenChannel', 1, self.warn)
            else
                self.warn = calllib('mstir', 'mNewMexPrinter');
                self.warn_case = 0;
                calllib('mstir', 'mOpenChannel', 1, self.warn)
            end
            if nargin > 2
                if ~ischar(errr)
                    error('Printer:wrong_arg', '??? Wrong errr in Printer');
                elseif strcmpi(errr, 'stdout') ~= 0
                    self.errr = calllib('mstir', 'mNewMexPrinter');
                    self.errr_case = 0;
                else
                    self.errr = calllib('mstir', 'mNewTextWriter', errr);
                    self.errr_case = 1;
                end
                calllib('mstir', 'mOpenChannel', 2, self.errr)
            else
                self.errr = calllib('mstir', 'mNewMexPrinter');
                self.errr_case = 0;
                calllib('mstir', 'mOpenChannel', 2, self.errr)
            end
        end
        function delete(self)
            calllib('mstir', 'mCloseChannel', 0, self.info)
            if self.info_case == 0
                calllib('mstir', 'mDeleteMexPrinter', self.info)
            else
                calllib('mstir', 'mDeleteTextWriter', self.info)
            end
            calllib('mstir', 'mCloseChannel', 1, self.warn)
            if self.warn_case == 0
                calllib('mstir', 'mDeleteMexPrinter', self.warn)
            else
                calllib('mstir', 'mDeleteTextWriter', self.warn)
            end
            calllib('mstir', 'mCloseChannel', 2, self.errr)
            if self.errr_case == 0
                calllib('mstir', 'mDeleteMexPrinter', self.errr)
            else
                calllib('mstir', 'mDeleteTextWriter', self.errr)
            end
        end
    end
end