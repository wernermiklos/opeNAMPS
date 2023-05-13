function out = ABMPS_apply_singlesite_op(obj, oppos, op, oprep, sites, symmetries, fermionicity)
%ABMPS_APPLY_SINGLESITE_OP Summary of this function goes here
%   Detailed explanation goes here
    if ~ABMPS_is_left(obj)
        error('MPS should be left canonical!!')
    end


    out = obj;
    if fermionicity < 0
        for pos = 1:(oppos-1)
            out.left_matrices{pos} = NTdot(out.left_matrices{pos},sites{pos}.operators('ph'),{'tau'},{'tau'},{{},{{'tau~','tau'}}});
        end
    end
    out.left_matrices{oppos} = NTdot(out.left_matrices{oppos},op,{'tau'},{'tau'},{{},{{'tau~','tau'}}});
    
    
    % From here we update the quantum numbers of the MPS matrices (this
    % could be done more beautifully, but it works for now. In the NA case
    % the Wigner-6J chain replaces this.
    [irreps, blocks, ~] = NTexport_all_data_blocks(out.left_matrices{oppos},{{'t_in',1},{'tau',1},{'t_out',1}},{'t_in','tau','t_out'});
    out.left_matrices{oppos} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},obj.no_of_symmetries);
    newirreps = cell(1,length(irreps));
    for bID = 1:length(irreps)
        newirreps{bID} = irreps{bID};
        tmp = SUPER_fusion_rule(symmetries,irreps{bID}{3},oprep);
        newirreps{bID}{3} = tmp{1}{1};
        out.left_matrices{oppos} = NTset_block(out.left_matrices{oppos},...
                                   {{'t_in',1},{'tau',1},{'t_out',1}},...
                                   newirreps{bID},...
                                   {'t_in','tau','t_out'},...
                                   blocks{bID});
    end
    [irreps, blocks, ~] = NTexport_all_data_blocks(out.schmidt_list{oppos},{{'t_left',1},{'t_right',1}},{'t_left', 't_right'});
    out.schmidt_list{oppos} = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},obj.no_of_symmetries);
    newirreps = cell(1,length(irreps));
    for bID = 1:length(irreps)
        newirreps{bID} = irreps{bID};
        tmp = SUPER_fusion_rule(symmetries,irreps{bID}{1},oprep);
        newirreps{bID}{1} = tmp{1}{1};
        newirreps{bID}{2} = SUPER_conj_rep(symmetries,tmp{1}{1});
        out.schmidt_list{oppos} = NTset_block(out.schmidt_list{oppos},...
                                              {{'t_left',1},{'t_right',1}},...
                                              newirreps{bID},...
                                              {'t_left','t_right'},...
                                              blocks{bID});
    end
    
    
    
    for pos = (oppos+1):(obj.chain_length)
        [irreps, blocks, ~] = NTexport_all_data_blocks(out.left_matrices{pos},{{'t_in',1},{'tau',1},{'t_out',1}},{'t_in','tau','t_out'});
        out.left_matrices{pos} = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},obj.no_of_symmetries);
        newirreps = cell(1,length(irreps));
        for bID = 1:length(irreps)
            newirreps{bID} = irreps{bID};
            tmp = SUPER_fusion_rule(symmetries,irreps{bID}{3},oprep);
            newirreps{bID}{3} = tmp{1}{1};
            tmp = SUPER_fusion_rule(symmetries,irreps{bID}{1},oprep);
            newirreps{bID}{1} = tmp{1}{1};
            out.left_matrices{pos} = NTset_block(out.left_matrices{pos},...
                                       {{'t_in',1},{'tau',1},{'t_out',1}},...
                                       newirreps{bID},...
                                       {'t_in','tau','t_out'},...
                                       blocks{bID});
        end
        [irreps, blocks, ~] = NTexport_all_data_blocks(out.schmidt_list{pos},{{'t_left',1},{'t_right',1}},{'t_left', 't_right'});
        out.schmidt_list{pos} = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},obj.no_of_symmetries);
        newirreps = cell(1,length(irreps));
        for bID = 1:length(irreps)
            newirreps{bID} = irreps{bID};
            tmp = SUPER_fusion_rule(symmetries,irreps{bID}{1},oprep);
            newirreps{bID}{1} = tmp{1}{1};
            newirreps{bID}{2} = SUPER_conj_rep(symmetries,tmp{1}{1});
            out.schmidt_list{pos} = NTset_block(out.schmidt_list{pos},...
                                                  {{'t_left',1},{'t_right',1}},...
                                                  newirreps{bID},...
                                                  {'t_left','t_right'},...
                                                  blocks{bID});
        end
        if pos == obj.chain_length
            out = ABMPS_set_schmidt(out,out.schmidt_list{pos});
        end
    end

