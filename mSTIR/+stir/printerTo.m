classdef printerTo
    properties
        channel
        printer
        type
    end
    methods
        function self = printerTo(dest, channel)
            if strcmpi(dest, 'stdout') ~= 0
                self.printer = calllib('mstir', 'mNewMexPrinter');
                self.type = 0;
            else
                self.printer = calllib('mstir', 'mNewTextWriter', dest);
                self.type = 1;
            end
            if nargin < 2
                calllib('mstir', 'mOpenChannel', 0, self.printer)
                calllib('mstir', 'mOpenChannel', 1, self.printer)
                calllib('mstir', 'mOpenChannel', 2, self.printer)
                self.channel = -1;
            else
                calllib('mstir', 'mOpenChannel', channel, self.printer)
                self.channel = channel;
            end
        end
        function delete(self)
            if self.type == 0
                calllib('mstir', 'mDeleteMexPrinter', self.printer)
            else
                calllib('mstir', 'mDeleteTextWriter', self.printer)
            end
            calllib('mstir', 'mCloseChannel', self.channel)
        end
    end
end