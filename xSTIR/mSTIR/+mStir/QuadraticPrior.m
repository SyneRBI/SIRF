classdef QuadraticPrior < mStir.Prior
    methods
        function self = QuadraticPrior()
            self.name = 'QuadraticPrior';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
%            calllib('mstir', 'mSTIR_deleteObject', self.handle)
            self.handle = [];
        end
    end
end