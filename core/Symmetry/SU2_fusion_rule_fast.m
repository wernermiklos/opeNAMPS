function out = SU2_fusion_rule_fast(dummy, Gamma1, Gamma2)
%SU2_FUSION_RULE_FAST Is used to replace the slow (Clebsch-based)
%fusion rule function handle in SU2sym.
Gamma_min = abs(Gamma1-Gamma2) + 1;
Gamma_max = Gamma1 + Gamma2 - 1;
out = cell(1,(Gamma_max - Gamma_min)/2 + 1);
i = 1;
for Gamma = Gamma_min:2:Gamma_max
    out{i} = {Gamma,1};
    i = i + 1;
end
end

