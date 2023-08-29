function out = ABSITE_U1_spinless_fermion(U1sym, Q_empty, Q_occupied)

%
%   Constructs Abelian site using the U1 symmetry for a spinless fermion. 
%
%   In order to be usable in infinite chain algorithms, the user can
%   specify the charge of the empty and occupied site freely by passing
%   integer values in Q_empty and Q_occupied.
%
%
%   Just for reminder, within the code we use the following representation
%   indices for the charges:
%   Representation | 1   2  3   4  5   6  7  ...
%          Charge Q| 0  -1  1  -2  2  -3  3
%             
%
%  Possible states
%
%   States    Q              R_Q
%     |0>     Q_empty        U1sym.Qnum_to_Gamma(Q_empty)     
%     |1>     Q_occupied     U1sym.Qnum_to_Gamma(Q_occupied)
%     
%

NO_OF_SYMMETRIES = 1;
rep_empty = U1sym.Qnum_to_Gamma(Q_empty);
rep_occupied = U1sym.Qnum_to_Gamma(Q_occupied);

fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{rep_empty,rep_occupied},{'tau~','tau'},[1]);

f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
f = NTset_block(f,{{'tau',1},{'tau~',1}},{rep_occupied,rep_empty},{'tau~','tau'},[1]);

ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{rep_empty,rep_empty},{'tau~','tau'},[1]);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{rep_occupied,rep_occupied},{'tau~','tau'},[-1]);

out = ABSITE_create(1,{{rep_empty,1},{rep_occupied,1}},'FERMIONS',1);


%Occupation number
n = NTdot(fdag, f, {'tau'},{'tau~'});

out.name = 'ABSITE_U1_spinless_fermion';

out = ABSITE_define_operator(out,'fdag',fdag,-1);
out = ABSITE_define_operator(out,'f',f,-1);
out = ABSITE_define_operator(out,'ph',ph,1);
out = ABSITE_define_operator(out,'n',n,1);

out = ABSITE_dot_define(out,'fdag','f','fdag_f');

end