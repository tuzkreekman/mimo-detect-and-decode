classdef sigmoidLayer < nnet.layer.Layer %nnet.cnn.layer.SoftmaxLayer
    properties 
    end

    methods
        function Z = predict(~,X)
            % Forward input data through the layer and output the result
            Z = exp(X)./(exp(X)+1);
        end
        function dLdX = backward(~,X, ~,dLdZ,~)
            % Backward propagate the derivative of the loss function through 
            % the layer 
            dLdX = (exp(X)./(exp(X)+1)).*(1- (exp(X)./(exp(X)+1))) .* dLdZ;
        end
    end
end


