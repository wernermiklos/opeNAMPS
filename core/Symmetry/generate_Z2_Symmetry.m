function Z2sym = generate_Z2_Symmetry()
%GENERATE_Z2_SYMMETRY Generates the Z2 (parity) symmetry
%   There are two representations: rep = 1:  trivial (even)
%                                  rep = 2:  nontrivial (odd)
CG_Z2 = NAtensor({'m1', 'm2', 'M', 'alpha'}, {'i', 'i', 'o', 'o'} , ...
        {[1], [2], [3], [1, 2, 3]}, 1);
CG_Z2 = NTset_block(CG_Z2,{{'m1',1},{'m2',1},{'M',1}},{[1],[1],[1]},{'m1','m2','M','alpha'}, [1]);
CG_Z2 = NTset_block(CG_Z2,{{'m1',1},{'m2',1},{'M',1}},{[1],[2],[2]},{'m1','m2','M','alpha'}, [1]);
CG_Z2 = NTset_block(CG_Z2,{{'m1',1},{'m2',1},{'M',1}},{[2],[1],[2]},{'m1','m2','M','alpha'}, [1]);
CG_Z2 = NTset_block(CG_Z2,{{'m1',1},{'m2',1},{'M',1}},{[2],[2],[1]},{'m1','m2','M','alpha'}, [1]);


Z2sym = SYM_create(CG_Z2,'Z2 Symmetry');

% We overwrite the selection rule and conj reps functions with their faster
% version.
Z2sym.selection_rule = @(dummy,Gamma1,Gamma2,Gamma) (mod(Gamma1 + Gamma2 + Gamma, 2) == 1);
Z2sym.conj_reps = @(Gamma) Gamma;
Z2sym.fusion_rule = @(dummy,Gamma1,Gamma2) {{mod(Gamma1 + Gamma2,2) + 1 ,1 }};


end

