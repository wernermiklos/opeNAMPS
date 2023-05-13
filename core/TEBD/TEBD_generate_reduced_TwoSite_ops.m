function TEBD_env = TEBD_generate_reduced_TwoSite_ops(TEBD_env)
%TEBD_GENERATE_TWOSITE_RED Summary of this function goes here
%   Detailed explanation goes here
    for pos = 1:length(TEBD_env.TwoSiteHlist_full)
        if pos ~= 1 && pos ~= length(TEBD_env.TwoSiteHlist_full) && abs(NTget_max_tensor_element(NTsubtr(TEBD_env.TwoSiteHlist_full{pos},TEBD_env.TwoSiteHlist_full{pos-1})))<10^(-10)
            TEBD_env.TwoSiteHlist_red{pos} = TEBD_env.TwoSiteHlist_red{pos-1};
        else
            TEBD_env.TwoSiteHlist_red{pos} = LNTdot(TEBD_env.envelope_list{pos},...
                                                    TEBD_env.TwoSiteHlist_full{pos}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'});
        end
        if pos ~= 1 && pos ~= length(TEBD_env.TwoSiteHlist_full) && abs(NTget_max_tensor_element(NTsubtr(TEBD_env.Ulist_full{pos},TEBD_env.Ulist_full{pos-1})))<10^(-10)
            TEBD_env.Ulist_red{pos} = TEBD_env.Ulist_red{pos-1};
        else
            TEBD_env.Ulist_red{pos} = LNTdot(TEBD_env.envelope_list{pos},...
                                                    TEBD_env.Ulist_full{pos}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'});
        end
        if pos ~= 1 && pos ~= length(TEBD_env.TwoSiteHlist_full) && abs(NTget_max_tensor_element(NTsubtr(TEBD_env.Uhalflist_full{pos},TEBD_env.Uhalflist_full{pos-1})))<10^(-10)
            TEBD_env.Uhalflist_red{pos} = TEBD_env.Uhalflist_red{pos-1};
        else
            TEBD_env.Uhalflist_red{pos} = LNTdot(TEBD_env.envelope_list{pos},...
                                                    TEBD_env.Uhalflist_full{pos}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'}, ...
                                                    {'mu_1','mu_2','mu_1~','mu_2~'});
        end
    end

end

