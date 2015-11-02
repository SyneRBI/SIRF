classdef printerTo
    properties
        channel_
        printer_
        type_
    end
    methods
        function self = printerTo(dest, channel)
            if strcmpi(dest, 'stdout') ~= 0
                self.printer_ = calllib('mstir', 'mNewMexPrinter');
                self.type_ = 0;
            else
                self.printer_ = calllib('mstir', 'mNewTextWriter', dest);
                self.type_ = 1;
            end
            if nargin < 2
                calllib('mstir', 'mOpenChannel', 0, self.printer_)
                calllib('mstir', 'mOpenChannel', 1, self.printer_)
                calllib('mstir', 'mOpenChannel', 2, self.printer_)
                self.channel_ = -1;
            else
                calllib('mstir', 'mOpenChannel', channel, self.printer_)
                self.channel_ = channel;
            end
        end
        function delete(self)
            if self.type_ == 0
                calllib('mstir', 'mDeleteMexPrinter', self.printer_)
            else
                calllib('mstir', 'mDeleteTextWriter', self.printer_)
            end
            calllib('mstir', 'mCloseChannel', self.channel_)
        end
    end
end