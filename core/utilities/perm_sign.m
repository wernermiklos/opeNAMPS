function out = perm_sign(perm)
%PERM_SIGN Summary of this function goes here
%   Returns the sign of permutation in perm. perm must contain a
%   permutation of numbers 1:N
A = eye(numel(perm));
out = det(A(:,perm));

end

