function out = ABMPS_create(chain_length, cut_position, no_of_symmetries)
%MPS_CREATE Creates an empty MPS object
%   chain_length:      chain length
%   cut_position:      the position of the last site in the "left"-part.
%   no_of_symmetries:  the number of symmetries
out.type = 'ABMPS';
out.left_matrices = cell(1,cut_position);
out.no_of_symmetries = no_of_symmetries;
out.cut_position = cut_position;
out.chain_length = chain_length;
out.schmidt_list = cell(1,chain_length);
for i = 1:length(out.left_matrices)
    out.left_matrices{i} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},out.no_of_symmetries);
end

out.right_matrices = cell(1,chain_length-cut_position);
for i = 1:length(out.right_matrices)
    out.right_matrices{i} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},out.no_of_symmetries);
end



out.schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},out.no_of_symmetries);


end

