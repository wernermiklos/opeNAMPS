function out = tensordot(A,B,contlegsA,contlegsB,freelegsA,freelegsB)
%TENSORDOT Summary of this function goes here
%   Detailed explanation goes here
    outshape = [size(A,freelegsA),size(B,freelegsB)];
    out = tensorprod(A,B,contlegsA,contlegsB);
    if length(outshape) >= 2
        out = reshape(out,outshape);
    end
end

