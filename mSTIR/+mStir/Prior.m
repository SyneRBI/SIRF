classdef Prior < handle
    properties
        name
        handle
    end
    methods
        function self = Prior()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_penalisation_factor(self, value)
            mStir.setParameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', value, 'f')
        end
        function value = get_penalisation_factor(self)
            value = mStir.parameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', 'f');
        end
    end
end