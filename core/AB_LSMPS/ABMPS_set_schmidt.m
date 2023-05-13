function out = ABMPS_set_schmidt(obj,schmidt)
%MPS_SET_SCHMIDT Summary of this function goes here
%   Detailed explanation goes here
out = obj;
out.schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},out.no_of_symmetries);  % we erase the original
out.schmidt = NTadd(out.schmidt,schmidt);
out.schmidt_list{out.cut_position} = out.schmidt;
end

