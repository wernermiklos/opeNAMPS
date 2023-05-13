function TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env, bond_pos, OpProd_list)
%ABTEBD_SET_TWOSITEH Summary of this function goes here
%   OpProd_list:  {{OpProd1, strength1}, {OpProd2, strength2}, ...}
%                   OpProd# should be generated using
%                   COUP_generate_TwoSite_full_OpProd
    H = 0;
    for i = 1:length(OpProd_list)
        H = NTadd(H,NTmult(OpProd_list{i}{1},OpProd_list{i}{2}));
    end

    TEBD_env.TwoSiteHlist_full{bond_pos} = H;





end

