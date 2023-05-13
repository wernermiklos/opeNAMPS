        
function out = ABSITE_spinhalf_fermion( varargin )
%
%   Constructs a general site based on the symmetries that are used
%
%
%   varargin can be a list with the names of the symmetries
%   Possible names to be used are 'U1_twoSz', 'U1_Q', 'Z1'
%   We can provide 0, 1 or 2 symmetries
%
%   Representation | 1   2  3   4  5   6  7  ...
%          Charge Q| 0  -1  1  -2  2  -3  3
%              2*Sz| 0  -1  1  -2  2  -3  3
%
%
%
%   States   2*Sz   Q   R_2Sz  R_Q
%     |0>   0     -1     1      2
%     |up>  1      0     3      1
%     |dn> -1      0     2      1
%  |up dn>  0      1     1      3
%
%
%
% . ↑ = char(8593)
% . ↓ = char(8595)

if nargin > 1
    error('ABSITE_spinhalf_fermion() must have 1 argument');
end

% Read the names of the symmetries that we plan to use and their number
if nargin == 1
    symmetry_names = varargin{1};
    no_of_symmetries =numel(symmetry_names);
end

%if no arguments are provided that consider the two symmetries case
if nargin == 0
    symmetry_names={'U1_twoSz', 'U1_Q'};
    no_of_symmetries =numel(symmetry_names);
end

%% Spin-half fermionic site with U(1)xU(1) symmetries for twoSz and Q

if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_twoSz')  && strcmp(symmetry_names{2},'U1_Q')
    
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
    
    
    out = ABSITE_create(2,{{[1,2],1},{[3,1],1}, {[2,1],1}, {[1,3],1}},'FERMIONS');
    
end

%% Spin-half fermionic site with U(1)xU(1) symmetries for Q and twoSz

if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q')  && strcmp(symmetry_names{2},'U1_twoSz')
    
    fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2,1],[1,3]},{'tau~','tau'},[1]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1,2],[3,1]},{'tau~','tau'},[1]);
   
    fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[2,1],[1,2]},{'tau~','tau'},[1]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1,3],[3,1]},{'tau~','tau'},[-1]);
    
    f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1,3],[2,1]},{'tau~','tau'},[1]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3,1],[1,2]},{'tau~','tau'},[1]);
   
    f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1,2],[2,1]},{'tau~','tau'},[1]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[3,1],[1,3]},{'tau~','tau'},[-1]);
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2,1],[2,1]},{'tau~','tau'},[1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,3],[1,3]},{'tau~','tau'},[-1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,2],[1,2]},{'tau~','tau'},[-1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3,1],[3,1]},{'tau~','tau'},[1]);
    
    
    out = ABSITE_create(2,{{[2,1],1},{[1,3],1}, {[1,2],1}, {[3,1],1}},'FERMIONS');
    
end


%% Spin-half fermionic site with U(1) symmetry for Q
if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'U1_Q')
    
    fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[1,0]');   
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[0,1]);    
 
    fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[0,1]');
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[-1,0]);
    
    f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[1,0]);   
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[0,1]');
    
    f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[0,1]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[-1,0]');
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2],[2]},{'tau~','tau'},[1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},-eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3],[3]},{'tau~','tau'},[1]);
    
    
    out = ABSITE_create(1,{{[2],1},{[1],2}, {[3],1}},'FERMIONS');
    
end

%% Spin-half fermionic site with U(1) symmetry for twoSz
if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'U1_twoSz')
    
    fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[0,1]');  
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[1,0]);
  
    fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[0,-1]');  
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[1,0]);
    
    f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[0,1]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[1,0]');
    
    f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[1,0]');
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[0,-1]);
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2],[2]},{'tau~','tau'},[-1]);   %single occ
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},+eye(2));   %double and empty
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3],[3]},{'tau~','tau'},[-1]);   %single occ
    
    
    out = ABSITE_create(1,{{[2],1},{[1],2}, {[3],1}},'FERMIONS');
    
end

%% Spin-half fermionic site with Z(1) symmetry  (corresponding to no symmetry at all)
if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'Z1')
    
    fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]);
    
    fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0 0 0 ; 0 0 0 0; 1 0 0 0; 0 0 -1 0]);
    
    
    f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]);
    
    f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0 1 0; 0 0 0 -1; 0 0 0 0; 0 0 0 0]);
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [1 0 0 0; 0 -1 0 0 ; 0 0 -1 0; 0 0 0 1]);
    
    
    out = ABSITE_create(1,{{[1],4}},'FERMIONS');
    
end

%% Construct the other operators and add them to the site.
% Spin-up occupation
n_u = NTdot(fdag_u, f_u, {'tau'},{'tau~'});

% Spin-down occupation
n_d = NTdot(fdag_d, f_d, {'tau'},{'tau~'});

%Total occupation
n = NTadd(NTdot(fdag_u, f_u, {'tau'},{'tau~'}), NTdot(fdag_d, f_d, {'tau'},{'tau~'}));

Sz = NTadd(NTmult(NTdot(fdag_u, f_u, {'tau'},{'tau~'}), 0.5), NTmult(NTdot(fdag_d, f_d, {'tau'},{'tau~'}),-0.5 ));
Sp = NTdot(fdag_u, f_d, {'tau'},{'tau~'});
Sm = NTdot(fdag_d, f_u, {'tau'},{'tau~'});


%Hubbard coupling: fdag_u fdag_d f_d f_u
hub = NTdot(fdag_u, NTdot(fdag_d, NTdot(f_d, f_u, {'tau'}, {'tau~'}), {'tau'}, {'tau~'}), {'tau'}, {'tau~'});


out.name = 'ABSITE_spinhalf_fermion';
out = ABSITE_define_operator(out,'fdag_u',fdag_u,-1);
out = ABSITE_define_operator(out,'fdag_d',fdag_d,-1);
out = ABSITE_define_operator(out,'f_u',f_u,-1);
out = ABSITE_define_operator(out,'f_d',f_d,-1);
out = ABSITE_define_operator(out,'ph',ph,1);
out = ABSITE_define_operator(out,'n_u',n_u,1);
out = ABSITE_define_operator(out,'n_d',n_d,1);
out = ABSITE_define_operator(out,'n',n,1);
out = ABSITE_define_operator(out,'hub',hub,1);

out = ABSITE_define_operator(out,'Sz',Sz,1);
out = ABSITE_define_operator(out,'Sp',Sp,1);
out = ABSITE_define_operator(out,'Sm',Sm,1);

end
