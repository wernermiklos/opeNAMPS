function out = NTsubtr(A, B)
% Function for subtraction: A - B. 
% Works only if the two tensors are consistent.
    out = NTadd(A, NTneg(B));
end

