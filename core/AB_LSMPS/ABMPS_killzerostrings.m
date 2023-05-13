function obj = ABMPS_killzerostrings(obj)
%ABMPS_KILLZEROS Kills blocks of MPS matrices that give zero contribution
%to the state (can be used in CAS initialization)
%   Detailed explanation goes here

    %left part
    t_out_irreps = NTget_irrep_set(obj.schmidt,{'t_left',1});
    for i = obj.cut_position:-1:1
        obj.left_matrices{i} = NTirrep_trunc_copy(obj.left_matrices{i},{{{'t_out',1},t_out_irreps}});
        if ~isempty(obj.schmidt_list{i})
            obj.schmidt_list{i} = NTirrep_trunc_copy(obj.schmidt_list{i},{{{'t_left',1},t_out_irreps}});
        end
        t_out_irreps = NTget_irrep_set(obj.left_matrices{i},{'t_in',1});
    end



    %right part
    t_out_irreps = NTget_irrep_set(obj.schmidt,{'t_right',1});
    for i = (obj.cut_position+1):obj.chain_length 
        obj.right_matrices{obj.chain_length-i+1} = NTirrep_trunc_copy(obj.right_matrices{obj.chain_length-i+1},{{{'t_out',1},t_out_irreps}});
        t_out_irreps = NTget_irrep_set(obj.right_matrices{obj.chain_length-i+1},{'t_in',1});
        if ~isempty(obj.schmidt_list{i})
            obj.schmidt_list{i} = NTirrep_trunc_copy(obj.schmidt_list{i},{{{'t_right',1},t_out_irreps}});
        end
    end
end

