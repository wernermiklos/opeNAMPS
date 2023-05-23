function out = MOD_ABSITE_gen_spinhalf_fermion_U1xU1()
%MOD_SITE_GEN_SPINLESS_FERMION_U1 Summary of this function goes here
%   Detailed explanation goes here
%   first U1:  spin-z
%   second U1: charge 

fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1,2],[3,1]},{'tau~','tau'},[1]);
fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2,1],[1,3]},{'tau~','tau'},[1]);
fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1,2],[2,1]},{'tau~','tau'},[1]);
fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[3,1],[1,3]},{'tau~','tau'},[-1]);
f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3,1],[1,2]},{'tau~','tau'},[1]);
f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1,3],[2,1]},{'tau~','tau'},[1]);
f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[2,1],[1,2]},{'tau~','tau'},[1]);
f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1,3],[3,1]},{'tau~','tau'},[-1]);


ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,2],[1,2]},{'tau~','tau'},[1]);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3,1],[3,1]},{'tau~','tau'},[-1]);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2,1],[2,1]},{'tau~','tau'},[-1]);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,3],[1,3]},{'tau~','tau'},[1]);

n = NTadd(NTdot(fdag_u, f_u, {'tau'},{'tau~'}), NTdot(fdag_d, f_d, {'tau'},{'tau~'}));

%Hubbard coupling: fdag_u fdag_d f_d f_u
hub = NTdot(fdag_u, NTdot(fdag_d, NTdot(f_d, f_u, {'tau'}, {'tau~'}), {'tau'}, {'tau~'}), {'tau'}, {'tau~'});



out = ABSITE_create(2,{{[1,2],1},{[3,1],1}, {[2,1],1}, {[1,3],1}},'FERMIONS',true);
out = ABSITE_define_operator(out,'fdag_u',fdag_u,-1);
out = ABSITE_define_operator(out,'fdag_d',fdag_d,-1);
out = ABSITE_define_operator(out,'f_u',f_u,-1);
out = ABSITE_define_operator(out,'f_d',f_d,-1);
out = ABSITE_define_operator(out,'ph',ph,1);
out = ABSITE_define_operator(out,'n',n,1);
out = ABSITE_define_operator(out,'hub',hub,1);

end