function out = LNTsubtr(obj, other)
            % Overloaded function of subtraction obj - other. 
            % Works only if the two tensors are consistent.
            out = LNTadd(obj,LNTneg(other));
end

