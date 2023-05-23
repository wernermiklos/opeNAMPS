function out = SITE_generate_SU2xU1_spinhalf_fermion(SU2sym,U1sym,varargin)
%SITE_GENERATE_S1_SITE_SU2 creates a spin half fermionic site
%   Detailed explanation goes here

 parameters = struct();
    % Default parameters
 parameters.CHARGE_OF_EMPTY_SITE = 0;
 parameters.CHARGE_PER_FERMION = +1;
    % Update according to varargin
 parameters = parameter_updater(parameters,varargin);
 
 empty_rep = [1,U1sym.Qnum_to_Gamma(parameters.CHARGE_OF_EMPTY_SITE)];
 single_rep = [2,U1sym.Qnum_to_Gamma(parameters.CHARGE_OF_EMPTY_SITE + parameters.CHARGE_PER_FERMION)];
 double_rep = [1,U1sym.Qnum_to_Gamma(parameters.CHARGE_OF_EMPTY_SITE + 2*parameters.CHARGE_PER_FERMION)];
 fdag_rep = [2,U1sym.Qnum_to_Gamma(parameters.CHARGE_PER_FERMION)];
 f_rep = [2,U1sym.Qnum_to_Gamma(-parameters.CHARGE_PER_FERMION)];
 triv_rep = [1,1];
 spin_rep = [3,1];   %S=1 (dim=3) Q=0 representation;

 out = SITE_create({SU2sym,U1sym},{{empty_rep,1},{single_rep,1},{double_rep,1}},'FERMIONS',true);

cdag_u = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{single_rep,empty_rep},{'mu~','mu','tau~','tau'},[1;0]);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{double_rep,single_rep},{'mu~','mu','tau~','tau'},[0,1]);
cdag_d = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{single_rep,empty_rep},{'mu~','mu','tau~','tau'},[0;1]);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{double_rep,single_rep},{'mu~','mu','tau~','tau'},[-1,0]);
c_u = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
c_u = NTset_block(c_u,{{'tau~',1},{'tau',1}},{empty_rep,single_rep},{'mu~','mu','tau~','tau'},[1,0]);
c_u = NTset_block(c_u,{{'tau~',1},{'tau',1}},{single_rep,double_rep},{'mu~','mu','tau~','tau'},[0;1]);
c_d = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
c_d = NTset_block(c_d,{{'tau~',1},{'tau',1}},{empty_rep,single_rep},{'mu~','mu','tau~','tau'},[0,1]);
c_d = NTset_block(c_d,{{'tau~',1},{'tau',1}},{single_rep,double_rep},{'mu~','mu','tau~','tau'},[-1;0]);

ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{empty_rep,empty_rep},{'mu~','mu','tau~','tau'},[1]);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{single_rep,single_rep},{'mu~','mu','tau~','tau'},-eye(2));
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{double_rep,double_rep},{'mu~','mu','tau~','tau'},[1]);

n = NTadd(NTdot(cdag_u, c_u, {'tau','mu'},{'tau~','mu~'}), NTdot(cdag_d, c_d, {'tau','mu'},{'tau~','mu~'}));

%Hubbard coupling: fdag_u fdag_d f_d f_u
hub = NTdot(cdag_u, NTdot(cdag_d, NTdot(c_d, c_u, {'tau','mu'}, {'tau~','mu~'}), ...
                          {'tau','mu'}, {'tau~','mu~'}), ...
            {'tau','mu'}, {'tau~','mu~'});

%Spin operators
% convention  S_i = 0.5* cdag \sigma_i c,  where cdag = (cdag_u cdag_d) and c = (c_u c_d), and sigma_i is the Pauli matrix (i = x,y,z) 
Sx = NTmult(NTadd(NTdot(cdag_u, c_d, {'tau','mu'},{'tau~','mu~'}), NTdot(cdag_d, c_u, {'tau','mu'},{'tau~','mu~'})),0.5);
Sy = NTmult(NTsubtr(NTdot(cdag_u, c_d, {'tau','mu'},{'tau~','mu~'}), NTdot(cdag_d, c_u, {'tau','mu'},{'tau~','mu~'})),-0.5i);
Sz = NTmult(NTsubtr(NTdot(cdag_u, c_u, {'tau','mu'},{'tau~','mu~'}), NTdot(cdag_d, c_d, {'tau','mu'},{'tau~','mu~'})),0.5);

S_1 = NTmult(NTadd(Sx,NTmult(Sy,1i)), -1/sqrt(2));
S_2 = Sz;
S_3 = NTmult(NTsubtr(Sx,NTmult(Sy,1i)), +1/sqrt(2));

out = SITE_define_tensor_operator(out,'fdag',{SU2sym,U1sym},fdag_rep,{cdag_u,cdag_d},'FERMIONICITY',-1);
out = SITE_define_tensor_operator(out,'f',{SU2sym,U1sym},f_rep,{c_d,NTneg(c_u)},'FERMIONICITY',-1); 
out = SITE_define_tensor_operator(out,'S',{SU2sym,U1sym},spin_rep,{S_1,S_2,S_3},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'ph',{SU2sym,U1sym},triv_rep,{ph},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'n',{SU2sym,U1sym},triv_rep,{n},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'hub',{SU2sym,U1sym},triv_rep,{hub},'FERMIONICITY',+1);
end
