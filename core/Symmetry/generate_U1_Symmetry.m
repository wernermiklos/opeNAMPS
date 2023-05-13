function U1sym = generate_U1_Symmetry(r1list, r2list, Rlist)
%GENERATE_U1_SYMMETRY Summary of this function goes here
%   Detailed explanation goes here
CG_U1 = generate_U1_CG(r1list, r2list, Rlist);
U1sym = SYM_create(CG_U1,'U1 Symmetry');

% We overwrite the selection rule and conj reps functions with their faster
% version.
U1sym.selection_rule = @(dummy,Gamma1,Gamma2,Gamma) ( floor(Gamma1/2)*(-1)^(Gamma1-1) +  floor(Gamma2/2)*(-1)^(Gamma2-1) ==  floor(Gamma/2)*(-1)^(Gamma-1));
U1sym.conj_reps = @(Gamma) max(1,Gamma + mod(Gamma+1,2)-mod(Gamma,2));
U1sym.fusion_rule = @(dummy,Gamma1,Gamma2) {{2*abs(floor(Gamma1/2)*(-1)^(Gamma1-1) +  floor(Gamma2/2)*(-1)^(Gamma2-1)) + ...
                                               (floor(Gamma1/2)*(-1)^(Gamma1-1) +  floor(Gamma2/2)*(-1)^(Gamma2-1) >= 0),1 }};
U1sym.Qnum_to_Gamma = @(Q) 2*abs(Q) + (Q >=0);
U1sym.Gamma_to_Qnum = @(Gamma) floor(Gamma/2)*(-1)^(Gamma-1);

end










