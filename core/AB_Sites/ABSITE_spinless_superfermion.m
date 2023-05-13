function out = ABSITE_spinless_superfermion( varargin )

%
%   Constructs a general site based on the charge symmetry
%
%   varargin{1}= symmetry_names
%   varargin{2} =convention for the creation/annihilation operators. 
%
%   varargin can be a list with the names of the symmetries
%   Possible names to be used are  'U1_Q', 'Z1'
%   We can provide 0, 1 or 2 symmetries
%
%   Representation | 1   2  3   4  5   6  7  ...
%          Charge Q| 0  -1  1  -2  2  -3  3
%
%
%  Possible states (4 in total)
%  The first state labes the regular one while the second
%  labels the tilde operator.
%
%
%   States             Q          R_Q
%----------------------------------------
%     |0,0>             0         1
%     |1,1>
%-----------------------------------------
%     |0,1>             -1        2
%-----------------------------------------
%     |1,0>              1        3
%-----------------------------------------


if nargin > 2
    error('ABSITE_spinless_fermion() must have at most 2 argument');
end

% Read the names of the symmetries that we plan to use and their number
if nargin==2
    symmetry_names = varargin{1};
    no_of_symmetries =numel(symmetry_names);
    convention = varargin{2};
end
if nargin == 1
    symmetry_names = varargin{1};
    no_of_symmetries =numel(symmetry_names);
    convention = 'real';
end

%if no arguments are provided that consider the two symmetries case
if nargin == 0
    symmetry_names={'U1_Q'};
    no_of_symmetries =numel(symmetry_names);
    convention ='real';
end


if  no_of_symmetries == 1 && strcmp(symmetry_names{1},'U1_Q')
    
    
    % regular convention for the creation/annihilation operators
    fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[1 0]);
    fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[0;1]);
    
    f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f = NTset_block(f,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[1;0]);
    f = NTset_block(f,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[0 1]);
    
    
    % one convention for the creation/annihilation operators.
    if strcmp(convention, 'real')
        
        tilde_fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[1 0]);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[0; -1]);
        
        tilde_f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_f =  NTset_block(tilde_f,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[0 -1]);
        tilde_f =  NTset_block(tilde_f,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[1;0]);
    end
    
    % Other possible convension for the tilde operators
    
    if strcmp(convention, 'complex')
        tilde_fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[1],[2]},{'tau~','tau'},[1.i 0]);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[3],[1]},{'tau~','tau'},[0; -1.i]);
        
        tilde_f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_f =  NTset_block(tilde_f,{{'tau',1},{'tau~',1}},{[1],[3]},{'tau~','tau'},[0 1.i]);
        tilde_f =  NTset_block(tilde_f,{{'tau',1},{'tau~',1}},{[2],[1]},{'tau~','tau'},[-1.i;0]);
        
    end
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2],[2]},{'tau~','tau'},[-1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3],[3]},{'tau~','tau'},[-1]);
    
    out = ABSITE_create(1,{{[1],2},{[2],1}, {[3], 1}},'FERMIONS');
    
end


if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'Z1')
    
    fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    fdag = NTset_block(fdag,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0 0 0 0; 0 0 0 1; 1 0 0 0; 0 0 0 0]);
    
    
    f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    f = NTset_block(f,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0 1 0; 0 0 0 0; 0 0 0 0; 0 1 0 0]);
    
    if strcmp(convention, 'real')
        tilde_fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0 0 0 0; 0 0 -1 0; 0 0 0 0; 1 0 0 0]);
        
        
        tilde_f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_f = NTset_block(f,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0 0 1; 0 0 0 0; 0 -1 0 0; 0 0 0 0]);
    end
    
    if strcmp(convention, 'complex')
        tilde_fdag = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_fdag = NTset_block(tilde_fdag,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0 0 0 0; 0 0 -1.i 0; 0 0 0 0; 1.i 0 0 0]);
        
        
        tilde_f = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
        tilde_f = NTset_block(f,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0 0 -1.i; 0 0 0 0; 0 1.i 0 0; 0 0 0 0]);
    end
    
    
    
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]);
    
    
    out = ABSITE_create(1,{{[1],4}},'FERMIONS');
    
end

%Occupation number
n = NTdot(fdag, f, {'tau'},{'tau~'});
tilde_n = NTdot(tilde_fdag, tilde_f, {'tau'},{'tau~'});

ddag = NTdot(fdag, tilde_fdag, {'tau'},{'tau~'});
d = NTdot(f, tilde_f, {'tau'},{'tau~'});

out.name = 'ABSITE_spinless_fermion';

out = ABSITE_define_operator(out,'fdag',fdag,-1);
out = ABSITE_define_operator(out,'f',f,-1);
out = ABSITE_define_operator(out,'n',n,1);

out = ABSITE_define_operator(out,'tilde_fdag',tilde_fdag,-1);
out = ABSITE_define_operator(out,'tilde_f',tilde_f,-1);
out = ABSITE_define_operator(out,'tilde_n',tilde_n,1);

out = ABSITE_define_operator(out,'ph',ph,1);
out = ABSITE_define_operator(out,'ddag',ddag,1);
out = ABSITE_define_operator(out,'d',d,1);


end