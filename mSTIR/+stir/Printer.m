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
            if isempty(info)
                self.info_case = -1;
            elseif strcmpi(info, 'stdout') ~= 0
                self.info = calllib('mstir', 'mNewMexPrinter');
                self.info_case = 0;
            else
                self.info = calllib('mstir', 'mNewTextWriter', info);
                self.info_case = 1;
            end
            if self.info_case ~= -1
                calllib('mstir', 'mOpenChannel', 0, self.info)
            else
                calllib('mstir', 'mCloseChannel', 0)
            end
%             self.warn = calllib('mstir', 'mNewMexPrinter');
%             self.warn_case = 0;
%             calllib('mstir', 'mOpenChannel', 1, self.warn)
%             self.errr = calllib('mstir', 'mNewMexPrinter');
%             self.errr_case = 0;
%             calllib('mstir', 'mOpenChannel', 2, self.errr)
%             if nargin < 2
%                 return
%             end
%             if isempty(warn)
%                 self.warn_case = -1;
%                 calllib('mstir', 'mCloseChannel', 1)
%             elseif strcmpi(warn, 'stdout') == 0
%                 self.warn = calllib('mstir', 'mNewTextWriter', warn);
%                 self.warn_case = 1;
%                 calllib('mstir', 'mOpenChannel', 1, self.warn)
%             end
%             if nargin < 3
%                 return
%             end
%             if isempty(errr)
%                 self.errr_case = -1;
%                 calllib('mstir', 'mCloseChannel', 2)
%             elseif strcmpi(errr, 'stdout') == 0
%                 self.errr = calllib('mstir', 'mNewTextWriter', errr);
%                 self.errr_case = 1;
%                 calllib('mstir', 'mOpenChannel', 2, self.errr)
%             end
        end
        function delete(self)
            if self.info_case ~= -1
                %calllib('mstir', 'mCloseChannel', 0)
                if self.info_case == 0
                    calllib('mstir', 'mDeleteMexPrinter', self.info)
                else
                    calllib('mstir', 'mDeleteTextWriter', self.info)
                end
            end
%             if self.warn_case ~= -1
%                 calllib('mstir', 'mCloseChannel', 1)
%                 if self.warn_case == 0
%                     calllib('mstir', 'mDeleteMexPrinter', self.warn)
%                 else
%                     calllib('mstir', 'mDeleteTextWriter', self.warn)
%                 end
%             end
%             if self.errr_case ~= -1
%                 calllib('mstir', 'mCloseChannel', 2)
%                 if self.errr_case == 0
%                     calllib('mstir', 'mDeleteMexPrinter', self.errr)
%                 else
%                     calllib('mstir', 'mDeleteTextWriter', self.errr)
%                 end
%             end
        end
    end
end