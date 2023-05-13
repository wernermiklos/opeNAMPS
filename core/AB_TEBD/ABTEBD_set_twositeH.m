function ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env, bond_pos, coupling_list)
%ABTEBD_SET_TWOSITEH Summary of this function goes here
%   coupling list:  {{'op1', 'op2', strength}, {'op1','op2', strength}, ...}
    H = 0;
    site1 = ABTEBD_env.Sites{bond_pos};
    site2 = ABTEBD_env.Sites{mod(bond_pos,ABTEBD_env.chain_length)+1};
    for i = 1:length(coupling_list)
        op1 = NTrename_legs(site1.operators(coupling_list{i}{1}),{'tau','tau~'},{'tau_1','tau_1~'});
        op2 = NTrename_legs(site2.operators(coupling_list{i}{2}),{'tau','tau~'},{'tau_2','tau_2~'});
        if ABTEBD_env.fermions
            if site1.fermionicities(coupling_list{i}{1}) * ...
                    site2.fermionicities(coupling_list{i}{2}) < 0
                error('the couplings must contain even number of fermi operators')
            end

            if site1.fermionicities(coupling_list{i}{1}) < 0
                op1 = NTdot(site1.operators('ph'),op1,{'tau~'},{'tau_1'},{{{'tau','tau_1'}},{}});
            end
        end
        H = NTadd(H, NTmult(NTsimple_mult(op1,op2), coupling_list{i}{3}));
    end
    ABTEBD_env.TwoSiteHlist{bond_pos} = H;


end

