function out = SITE_generate_U1_spinless_fermion(U1sym, Q_empty, Q_occupied)
%SITE_GENERATE_S1_SITE_SU2 Constructs a U1 spinless fermion site but using the general (non-Abelian) structures. 
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
     rep_fdag = U1sym.Qnum_to_Gamma(Q_occupied-Q_empty);       % representation of the creation operator.
     rep_f = U1sym.Qnum_to_Gamma(Q_empty-Q_occupied);          % representation of the annihilation operator.
     rep_triv = 1;
      
     
     out = SITE_create({U1sym},{{rep_empty,1},{rep_occupied,1}},'FERMIONS',true);
     
     fdag = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);
     fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{rep_empty,rep_occupied},{'tau~','mu~','tau','mu'},[1]);
        
     f =  NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);
     f = NTset_block(f,{{'tau',1},{'tau~',1}},{rep_occupied,rep_empty},{'tau~','mu~','tau','mu'},[1]);
        
     ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);
     ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{rep_empty,rep_empty},{'tau~','mu~','tau','mu'},[1]);
     ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{rep_occupied,rep_occupied},{'tau~','mu~','tau','mu'},[-1]);
    
     n = NTdot(fdag, f, {'tau','mu'},{'tau~','mu~'});
    
     out = SITE_define_tensor_operator(out,'fdag',{U1sym},rep_fdag,{fdag},'FERMIONICITY',-1);
     out = SITE_define_tensor_operator(out,'f',{U1sym},rep_f,{f},'FERMIONICITY',-1); 
     out = SITE_define_tensor_operator(out,'ph',{U1sym},rep_triv,{ph},'FERMIONICITY',+1);
     out = SITE_define_tensor_operator(out,'n',{U1sym},rep_triv,{n},'FERMIONICITY',+1);
end
