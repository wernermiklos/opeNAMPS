function Z1sym = generate_dummy_Z1_Symmetry()
%GENERATE_U1_SYMMETRY Summary of this function goes here
%   Detailed explanation goes here
CG_Z1 = generate_U1_CG([1], [1], [1]);
Z1sym = SYM_create(CG_Z1,'Z1 Symmetry');

% We overwrite the selection rule and conj reps functions with their faster
% version.
Z1sym.selection_rule = @(dummy,Gamma1,Gamma2,Gamma) (Gamma1 == 1) && (Gamma2 == 1) && (Gamma == 1);
Z1sym.conj_reps = @(Gamma) Gamma;
Z1sym.fusion_rule = @(dummy,Gamma1,Gamma2) {{1,(Gamma1 == 1) && (Gamma2 == 1)}};
Z1sym.Qnum_to_Gamma = @(Q) Q;
Z1sym.Gamma_to_Qnum = @(Gamma) Gamma;

end










