function outstruct = ABTEBD_svd_OLD(psi, symmetries, Mmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    outstruct = struct();
    outstruct.normsq_untr = NTdot(psi,NTconj(psi), ...
                                  {'t_left','tau_1','tau_2','t_right'}, ...
                                  {'t_left','tau_1','tau_2','t_right'});
    svdten = NAtensor({'t_left','tau_1','tau_2', 't_right','aux1','aux2'},{'i','i','i','i','o','i'},{[1],[2],[3],[4],[5],[5]},psi.no_of_symmetries);

    [irreps, blocks] = NTexport_all_data_blocks(psi,{{'t_left',1},{'tau_1',1},{'tau_2',1},{'t_right',1}},{'t_left','tau_1','tau_2', 't_right'});
    for bID = 1:length(blocks)
        Gamma_left = irreps{bID}{1};
        Gamma_1 = irreps{bID}{2};
        fusion = SUPER_fusion_rule(symmetries,Gamma_left,Gamma_1);
        Gamma_center = fusion{1}(1);
        svdten = NTset_block(svdten,...
                             {{'t_left',1},{'tau_1',1},{'tau_2',1},{'t_right',1},{'aux1',1}}, ...
                              [irreps{bID},Gamma_center], ...
                              {'t_left','tau_1','tau_2', 't_right','aux1','aux2'}, ...
                             blocks{bID});
    end
    [U,newschmidt,V,normsq] = NTsvd(svdten,{'t_left', 'tau_1', 'aux1'},{'aux2', 'tau_2', 't_right'},{'t_bond_left', 't_bond_right~'},{'i','o'},Mmax);
    conj_switch_tr = NAtensor({'t_bond_right','t_bond_right~'},{'i','i'},{[1],[2]},svdten.no_of_symmetries);
    bond_sectors = NTget_leg_sectors(newschmidt,'t_bond_left');
    for secID = 1:length(bond_sectors)
        Gamma_left = bond_sectors{secID}{1}{1};
        Gamma_right = SUPER_conj_rep(symmetries,Gamma_left);
        secdim = bond_sectors{secID}{2};
        conj_switch_tr = NTset_block(conj_switch_tr,{{'t_bond_right~',1},{'t_bond_right',1}},{Gamma_left, Gamma_right},{'t_bond_right~','t_bond_right'},eye(secdim));
    end
    outstruct.schmidt = NTdot(newschmidt, conj_switch_tr, {'t_bond_right~'},{'t_bond_right~'},{{{'t_bond_left','t_left'}},{{'t_bond_right','t_right'}}});
    outstruct.A = NTrename_legs(NTcombine_legs(U,{'aux1','t_bond_left'},'t_out'),{'t_left','tau_1'},{'t_in','tau'});
    outstruct.B = NTdot(NTcombine_legs(V, {'aux2', 't_bond_right~'},'t_bond_right~'), NTconj(conj_switch_tr), {'t_bond_right~'},{'t_bond_right~'}, ... 
              {{{'tau_2','tau'},{'t_right','t_in'}},{{'t_bond_right','t_out'}}});
    outstruct.normsq = normsq;
          


end

