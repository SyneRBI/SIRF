classdef MessageRedirector < handle
    properties
        info
        info_case
        warn
        warn_case
        errr
        errr_case
    end
    methods
        function self = MessageRedirector(info, warn, errr)
%         MessageRedirector(info, warn, errr) re-defines printing destinations.
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
            if ~ischar(info)
                error('MessageRedirector:wrong_arg', ...
                    '??? Wrong info in MessageRedirector consttructor');
            elseif strcmpi(info, 'stdout') ~= 0
                self.info = calllib('mstir', 'mNewMexPrinter');
                self.info_case = 0;
            else
                self.info = calllib('mstir', 'mNewTextWriter', info);
                self.info_case = 1;
            end
            calllib('mstir', 'mOpenChannel', 0, self.info)
            if nargin > 1
                if ~ischar(warn)
                    error('MessageRedirector:wrong_arg', ...
                        '??? Wrong warn in MessageRedirector constructor');
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
                    error('MessageRedirector:wrong_arg', ...
                        '??? Wrong errr in MessageRedirector constructor');
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
                h = calllib('mstir', 'mDeleteMexPrinter', self.info);
            else
                h = calllib('mstir', 'mDeleteTextWriter', self.info);
            end
            mUtilities.check_status('MessageRedirector:delete', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
            calllib('mstir', 'mCloseChannel', 1, self.warn)
            if self.warn_case == 0
                h = calllib('mstir', 'mDeleteMexPrinter', self.warn);
            else
                h = calllib('mstir', 'mDeleteTextWriter', self.warn);
            end
            mUtilities.check_status('MessageRedirector:delete', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
            calllib('mstir', 'mCloseChannel', 2, self.errr)
            if self.errr_case == 0
                h = calllib('mstir', 'mDeleteMexPrinter', self.errr);
            else
                h = calllib('mstir', 'mDeleteTextWriter', self.errr);
            end
            mUtilities.check_status('MessageRedirector:delete', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
    end
end