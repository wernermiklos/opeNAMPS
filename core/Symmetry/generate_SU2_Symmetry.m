function SU2sym = generate_SU2_Symmetry(r1list, r2list, Rlist, varargin)
%GENERATE_SU2_SYMMETRY Summary of this function goes here
%   Detailed explanation goes here

    parameters = struct();
    parameters.INFO = false;
    parameters = parameter_updater(parameters,varargin);

    CG_SU2 = generate_SU2_CG(r1list,r2list,Rlist,'INFO',parameters.INFO);

    SU2sym = SYM_create(CG_SU2,'SU2 symmetry');
    
    % Faster, general selection rule function.
    SU2sym.selection_rule = @(dummy,Gamma1,Gamma2,Gamma) (Gamma <= (Gamma1 + Gamma2 - 1)) &&...
                                                         (Gamma >= (abs(Gamma1-Gamma2) + 1)) &&...
                                                         (mod(Gamma,2) == mod(Gamma1 + Gamma2 - 1,2));
    SU2sym.conj_reps = @(Gamma) Gamma;
    SU2sym.fusion_rule = @(dummy,Gamma1,Gamma2) SU2_fusion_rule_fast(dummy,Gamma1,Gamma2);
    
    SU2sym.Qnum_to_Gamma = @(S) 2*S+1;
    SU2sym.Gamma_to_Qnum = @(Gamma) (Gamma-1)/2;

end

