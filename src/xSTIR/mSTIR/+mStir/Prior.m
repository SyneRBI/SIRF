classdef Prior < handle
    % Class for objects handling the prior: a penalty term to be added to the
    % objective function maximized by iterative reconstruction algorithms.
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
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_penalisation_factor(self, value)
%         Sets the factor by which the penalty term (prior) is to be multiplied
%         before adding to the objective function.
            mStir.setParameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', value, 'f')
        end
        function value = get_penalisation_factor(self)
%         Returns the penalty factor in front of the prior.
            value = mStir.parameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', 'f');
        end
        function grad = get_gradient(self, image)
%         Returns the value of the gradient of the prior for a given value of
%         the image.
            grad = mStir.ImageData();
            grad.handle = calllib('mstir', 'mSTIR_priorGradient', ...
                self.handle, image.handle);
            mUtil.checkExecutionStatus('ImageData', grad.handle)
        end
    end
end