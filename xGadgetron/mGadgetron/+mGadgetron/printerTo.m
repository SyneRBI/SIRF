classdef printerTo < handle
    properties
        channel
        printer
        type
    end
    methods
        function self = printerTo(dest, channel)
            lib = 'mutilities';
            if strcmpi(dest, 'stdout') ~= 0
                self.printer = calllib(lib, 'mNewMexPrinter');
                self.type = 0;
            else
                self.printer = calllib(lib, 'mNewTextWriter', dest);
                self.type = 1;
            end
            if nargin < 2
                calllib(lib, 'mOpenChannel', 0, self.printer)
                calllib(lib, 'mOpenChannel', 1, self.printer)
                calllib(lib, 'mOpenChannel', 2, self.printer)
                self.channel = -1;
            else
                calllib(lib, 'mOpenChannel', channel, self.printer)
                self.channel = channel;
            end
            %fprintf('text printer to channel %d of type %d created\n', self.channel, self.type)
        end
        function print(self, text)
            calllib('mutilities', 'mPrintText', text)
        end
        function delete(self)
            lib = 'mutilities';
            %fprintf('deleting text printer of type %d...', self.type)
            if self.type == 0
                calllib(lib, 'mDeleteMexPrinter', self.printer)
            else
                calllib(lib, 'mDeleteTextWriter', self.printer)
            end
            calllib(lib, 'mCloseChannel', self.channel, self.printer)
            %fprintf('ok\n')
        end
    end
end