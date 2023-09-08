function out = ABSITE_spinless_fermion( varargin )

%
%   Constructs a general site based on the charge symmetry 
%
%
%   varargin can be a list with the names of the symmetries
%   Possible names to be used are  'U1_Q', 'Z1'
%   We can provide 0, 1 or 2 symmetries
%
%   Representation | 1   2  3   4  5   6  7  ...
%          Charge Q| 0  -1  1  -2  2  -3  3
%             
%
%  Possible states
%
%   States    Q    R_Q
%     |0>     0     1      
%     |1>     1     3     
%     
%


if nargin > 1
     error('ABSITE_spinless_fermion() must have 1 argument');
end

% Read the names of the symmetries that we plan to use and their number
if nargin == 1
    symmetry_names = varargin{1};
    no_of_symmetries =numel(symmetry_names);
end

%if no arguments are provided that consider the two symmetries case
if nargin == 0
    symmetry_names={'U1_Q'};
    no_of_symmetries =numel(symmetry_names);
end

if  no_of_symmetries == 1 && strcmp(symmetry_names{1},'U1_Q')

fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[1]);

f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
f = NTset_block(f,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[1]);

ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},[1]);
ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3],[3]},{'tau~','tau'},[-1]);

out = ABSITE_create(1,{{[1],1},{[3],1}},'FERMIONS',1);

end


if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'Z1')
    
    fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0 0; 1 0]);
    
    
    f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f = NTset_block(f,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 1; 0 0]);
    
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [1 0 ; 0 -1]);
    
    
    out = ABSITE_create(1,{{[1],2}},'FERMIONS');
    
end

%Occupation number
n = NTdot(fdag, f, {'tau'},{'tau~'});

out.name = 'ABSITE_spinless_fermion';

out = ABSITE_define_operator(out,'fdag',fdag,-1);
out = ABSITE_define_operator(out,'f',f,-1);
out = ABSITE_define_operator(out,'ph',ph,1);
out = ABSITE_define_operator(out,'n',n,1);

out = ABSITE_dot_define(out,'fdag','f','fdag_f');

end