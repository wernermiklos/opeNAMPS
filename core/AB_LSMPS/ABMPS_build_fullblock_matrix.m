function outmatrix = ABMPS_build_fullblock_matrix(symmetries, t_in_sectors, tau_sectors, keptirreps_tau, keptirreps_t_out)
%ABMPS_BUILD_FULLBLOCK_MATRIX Creates an ABELIAN MPS matrix where all states are kept
%   t_in_sectors:       sectors of the the input block (without the site).
%                       {{{irreps},dim},...}
%   tau_sectors:        the new site
%   keptirreps_t_out:   we can choose which irreps are kept in the bond. If not given, all
%                       irreps are kept.
%   keptirreps_tau      we can choose which irreps are kept in the site. If not given, all
%                       irreps are kept.

    %t_in_map = containers.Map();
    t_out_map = containers.Map();
    %tau_map = containers.Map();
    
    t_out_tmp = containers.Map();

    blocks = {};
    for secID_t_in = 1:length(t_in_sectors)
        Gamma_in = t_in_sectors{secID_t_in}{1}{1};
        t_in_dim = t_in_sectors{secID_t_in}{2};
        for secID_tau = 1:length(tau_sectors)
            Gamma_loc = tau_sectors{secID_tau}{1}{1};
            kept = true;
            if ~isempty(keptirreps_tau)
                if any(cellfun(@(x) isequal(x,Gamma_loc),keptirreps_tau))
                    tau_dim = tau_sectors{secID_tau}{2};
                    %tau_map(char(Gamma_loc)) = tau_dim;
                else
                    kept = false;
                end
            else
                tau_dim = tau_sectors{secID_tau}{2};
                %tau_map(char(Gamma_loc)) = tau_dim;
            end
            if kept
                 outsecs = SUPER_fusion_rule(symmetries,Gamma_in,Gamma_loc);
                 for secID_out = 1:length(outsecs)
                     Gamma_out = outsecs{secID_out}{1};
                     kept_out = true;
                     if ~isempty(keptirreps_t_out)
                            if ~any(cellfun(@(x) isequal(x,Gamma_out),keptirreps_t_out))
                                kept_out = false;
                            end
                     end
                     if kept_out
                         blocks{end+1} = {{Gamma_in, Gamma_loc, Gamma_out}, t_in_dim, tau_dim};
                         if ~t_out_tmp.isKey(char(Gamma_out))
                             t_out_map(char(Gamma_out)) = t_in_dim*tau_dim;
                             t_out_tmp(char(Gamma_out)) = 0;
                         else
                             t_out_map(char(Gamma_out)) = t_out_map(char(Gamma_out)) + t_in_dim*tau_dim;
                         end
                     end
                    
                 end
            end
        end
    end

    outmatrix = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},length(symmetries));
    
    for bID = 1:length(blocks)
        Gamma_in = blocks{bID}{1}{1};
        Gamma_loc = blocks{bID}{1}{2};
        Gamma_out = blocks{bID}{1}{3};
        t_in_dim = blocks{bID}{2};
        tau_dim = blocks{bID}{3};
        new_block = zeros(t_in_dim, tau_dim, t_out_map(char(Gamma_out)));
        for t_in=1:t_in_dim
            for tau = 1:tau_dim
                    t_out_tmp(char(Gamma_out)) = t_out_tmp(char(Gamma_out)) + 1;
                    new_block(t_in,tau,t_out_tmp(char(Gamma_out))) = 1;
            end
        end
        outmatrix = NTset_block(outmatrix, ...
                                    {{'t_in',1},{'tau',1},{'t_out',1}}, ...
                                    {Gamma_in,Gamma_loc,Gamma_out}, ...
                                    {'t_in','tau','t_out'}, ...
                                    new_block);


    end
    


end

