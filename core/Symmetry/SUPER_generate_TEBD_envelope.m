function out = SUPER_generate_TEBD_envelope(Symmetries,Site1,Site2,varargin)
%SUPER_GENERATE_TEBD_ENVELOPE Summary of this function goes here
%   Detailed explanation goes here
parameters = struct();
parameters.('BONDREPS_PER_SYM') = {};
parameters = parameter_updater(parameters,varargin);

site1_reps_1sym = NTget_irrep_set_1sym(Site1.operators('id'),{'tau',1});
site2_reps_1sym = NTget_irrep_set_1sym(Site2.operators('id'),{'tau',1});
envelope_list = cell(1,length(Symmetries));
for symID = 1:length(Symmetries)
    CGtrunc1 = NTirrep_trunc_copy(Symmetries{symID}.CGtensor,{{{'m2',1},site1_reps_1sym{symID}}});
    CGtrunc2 = NTirrep_trunc_copy(Symmetries{symID}.CGtensor,{{{'m2',1},site2_reps_1sym{symID}}});
    if ~isempty(parameters.BONDREPS_PER_SYM)
        if ~isempty(parameters.BONDREPS_PER_SYM{symID})
            CGtrunc1 = NTirrep_trunc_copy(CGtrunc1,{{{'m1',1},parameters.BONDREPS_PER_SYM{symID}}, ...
                                                    {{'M',1}, parameters.BONDREPS_PER_SYM{symID}}});
            CGtrunc2 = NTirrep_trunc_copy(CGtrunc2,{{{'m1',1},parameters.BONDREPS_PER_SYM{symID}}, ...
                                                    {{'M',1}, parameters.BONDREPS_PER_SYM{symID}}});
        end
    end
    tmp = NTdot(CGtrunc1,CGtrunc2,{'M'},{'m1'},{{{'m1','m_left'},{'m2','mu_1'},{'alpha','alpha_1'}},...
                                                {{'m2','mu_2'},{'alpha','alpha_2'}}});
    tmp2 = NTdot(tmp,Symmetries{symID}.one_per_dim,{'M'},{'M~'});
    tmp = NTprime_all_names(NTconj(tmp));
    envelope_list{symID} = NTdot(tmp2,tmp,{'m_left','M'},{'m_left~','M~'});
end

out = envelope_list{1};
for symID = 2:length(Symmetries)
    out = LNTirrep_kron(out,envelope_list{symID});
end



end

