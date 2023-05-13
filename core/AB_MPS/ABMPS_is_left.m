function out = ABMPS_is_left(obj)
%ABMPS_IS_LEFT Summary of this function goes here
%   Returns true if ABMPS is fully left canonical
if isempty(obj.right_matrices)
    out = true;
else
    out = false;

end

