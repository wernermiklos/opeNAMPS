function out = ABMPS_set_matrix(obj,pos,matrix)
%MPS_SET_MATRIX Summary of this function goes here
%   Detailed explanation goes here
out = obj;
if pos <= obj.cut_position
    i = pos;
    % first erase
    out.left_matrices{i} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},out.no_of_symmetries);
    % then fill
    out.left_matrices{i} = NTadd(out.left_matrices{i},matrix);
else
    i = obj.chain_length-pos+1;
    % first erase
    out.right_matrices{i} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},out.no_of_symmetries);
    % then fill
    out.right_matrices{i} = NTadd(out.right_matrices{i},matrix);
end

end

