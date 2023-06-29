function copy = ABMPS_copy(obj)
%ABMPS_COPY Summary of this function goes here
%   Detailed explanation goes here
copy=ABMPS_create(obj.chain_length,obj.cut_position,obj.no_of_symmetries);

copy.left_matrices=obj.left_matrices;
copy.right_matrices=obj.right_matrices;
copy.schmidt_list=obj.schmidt_list;
copy.schmidt=obj.schmidt;

end

