classdef QuadraticPrior < stir.GeneralisedPrior
    methods
        function self = QuadraticPrior()
            self.name = 'QuadraticPrior';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mstir', 'mSTIR_deleteObject', self.handle, 'Prior')
        end
    end
end