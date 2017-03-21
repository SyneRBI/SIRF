classdef QuadraticPrior < mStir.Prior
    % Class for the prior that is a quadratic function of the image values.
    methods
        function self = QuadraticPrior()
            self.name = 'QuadraticPrior';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
            self.handle = [];
        end
    end
end