function out = ABSITE_gen_fermion_site(symmetries, rep)
%MOD_SITE_GEN_SPINLESS_FERMION_U1 Defines a generic fermion site (2
%dimensional)
%   symmetries: list of symmetries used (needed to define "operator_irreps"
%               properly.
%   rep:        representation (U(1) charges & Z_2 parity) of the filled
%               state. Empty state is always defined in the trivial representation.

    trivrep = ones(1,length(rep));
    conjrep = SUPER_conj_rep(symmetries,rep);
    no_of_symmetries = length(rep);
    fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},no_of_symmetries);
    fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{trivrep,rep},{'tau~','tau'},[1]);
    f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},no_of_symmetries);
    f = NTset_block(f,{{'tau',1},{'tau~',1}},{rep,trivrep},{'tau~','tau'},[1]);

    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},no_of_symmetries);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{trivrep,trivrep},{'tau~','tau'},[1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{rep,rep},{'tau~','tau'},[-1]);

    out = ABSITE_create(no_of_symmetries,{{trivrep,1},{rep,1}},'FERMIONS',true);
    out = ABSITE_define_operator_with_oprep(out,'fdag',fdag, rep, 'FERMIONICITY', -1, 'SYMMETRIES',symmetries);
    out = ABSITE_define_operator_with_oprep(out,'f', f, conjrep, 'FERMIONICITY', -1, 'SYMMETRIES',symmetries);
    out = ABSITE_define_operator_with_oprep(out,'ph',ph, trivrep, 'FERMIONICITY', 1, 'SYMMETRIES',symmetries);
    out = ABSITE_dot_define(out,'fdag','f','fdag_f','SYMMETRIES',symmetries);  
end