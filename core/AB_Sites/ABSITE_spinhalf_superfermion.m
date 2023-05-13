        
function out = ABSITE_spinhalf_superfermion( varargin )
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
%   States             Q      2*Sz      R_Q   R_2Sz
%----------------------------------------------------
%     |0,(↑↓) >        -2      0         4       1
%----------------------------------------------------     
%     |0,↑>             -1    -1         2       2
%     |↓, (↑↓)>         
%----------------------------------------------------  
%     |0, ↓>            -1     1         2       3
%     |↑, (↑↓)>         
%----------------------------------------------------
%     |0, 0>
%     |(↑↓), (↑↓)>      0      0         1       1
%     |↑, ↑>
%     |↓,↓ >
%----------------------------------------------------
%     |↓, 0>            1       -1       3      2
%     |(↑↓),↑>
%---------------------------------------------------
%     |↑,↓>             0       2         1     5
%----------------------------------------------------
%     |↓,↑>             0       -2        1     4
%---------------------------------------------------
%     |↑,0 >            1       1         3     3
%     |(↑↓),↓ >         
%---------------------------------------------------
%     |(↑↓), 0>         2       0         5     1
%---------------------------------------------------
% . ↑ = char(8593)
% . ↓ = char(8595)
%


if nargin > 2
    error('ABSITE_spinhalf_fermion() must have at most 2 argument');
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
    symmetry_names={'U1_Q', 'U1_twoSz'};
    no_of_symmetries =numel(symmetry_names);
    convention ='real';
end




%% Spin-half fermionic site with U(1)xU(1) symmetries for Q and twoSz

if   no_of_symmetries == 2 && strcmp(symmetry_names{1}, 'U1_Q')  && strcmp(symmetry_names{2},'U1_twoSz') && strcmp(convention, 'real')
    
    fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[4,1],[2,3]},{'tau~','tau'},[0; 1]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2,2],[1,1]},{'tau~','tau'},[0 0; 0 1; 1 0; 0 0]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[2,3],[1,5]},{'tau~','tau'},[1 0]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1,1],[3,3]},{'tau~','tau'},[1 0 0 0; 0 0 0 1 ]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[3,2],[5,1]},{'tau~','tau'},[ 1 0]);
    fdag_u = NTset_block(fdag_u,{{'tau',1},{'tau~',1}},{[1,4],[3,2]},{'tau~','tau'},[0; 1]);
     
    
    fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[4,1],[2,2]},{'tau~','tau'},[0; 1]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[2,2],[1,4]},{'tau~','tau'},[1 0]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[2,3],[1,1]},{'tau~','tau'},[0 0; 0 -1; 0 0; 1 0]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1,1],[3,2]},{'tau~','tau'},[1 0 0 0; 0 0 -1 0]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[3,3],[5,1]},{'tau~','tau'},[-1 0]);
    fdag_d = NTset_block(fdag_d,{{'tau',1},{'tau~',1}},{[1,5],[3,3]},{'tau~','tau'},[0;-1]);
    
   
    tilde_fdag_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[2,3],[4,1]},{'tau~','tau'},[1 0]);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[1,1],[2,2]},{'tau~','tau'},[1 0 0  0; 0 0 0 -1]);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[3,2],[1,4]},{'tau~','tau'},[-1 0]);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[3,3],[1,1]},{'tau~','tau'},[0 0; 0 1; -1 0; 0 0 ]);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[5,1],[3,2]},{'tau~','tau'},[0; 1]);
    tilde_fdag_u = NTset_block(tilde_fdag_u,{{'tau',1},{'tau~',1}},{[1,5],[2,3]},{'tau~','tau'},[0; -1]);
    
    
    tilde_fdag_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[2,2],[4,1]},{'tau~','tau'},[-1 0]);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[1,1],[2,3]},{'tau~','tau'},[1 0 0 0; 0 0 1 0]);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[3,2],[1,1]},{'tau~','tau'},[0 0; 0 -1; 0 0; -1 0 ]);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[3,3],[1,5]},{'tau~','tau'},[-1 0]);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[5,1],[3,3]},{'tau~','tau'},[0; 1]);
    tilde_fdag_d = NTset_block(tilde_fdag_d,{{'tau',1},{'tau~',1}},{[1,4],[2,2]},{'tau~','tau'},[0; 1]);
    
    
    f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[2,3],[4,1]},{'tau~','tau'},[0 1]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1,1],[2,2]},{'tau~','tau'},[0 0 1 0; 0 1 0 0]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3,2],[1,4]},{'tau~','tau'},[0 1]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[3,3],[1,1]},{'tau~','tau'},[1 0; 0 0; 0 0; 0 1]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[5,1],[3,2]},{'tau~','tau'},[1; 0]);
    f_u = NTset_block(f_u,{{'tau',1},{'tau~',1}},{[1,5],[2,3]},{'tau~','tau'},[1; 0]); 
       
    
    f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[2,2],[4,1]},{'tau~','tau'},[0 1]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1,1],[2,3]},{'tau~','tau'},[0 0 0 1; 0 -1 0 0]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[3,2],[1,1]},{'tau~','tau'},[1 0; 0 0; 0 -1; 0 0]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[3,3],[1,5]},{'tau~','tau'},[0 -1]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[5,1],[3,3]},{'tau~','tau'},[-1; 0]);
    f_d = NTset_block(f_d,{{'tau',1},{'tau~',1}},{[1,4],[2,2]},{'tau~','tau'},[1; 0]);
    
    
    tilde_f_u = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[4,1],[2,3]},{'tau~','tau'},[1; 0]);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[2,2],[1,1]},{'tau~','tau'},[1 0; 0 0; 0 0; 0 -1]);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[2,3],[1,5]},{'tau~','tau'},[0 -1]);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[1,1],[3,3]},{'tau~','tau'},[ 0 0 -1 0; 0 1 0 0]);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[3,2],[5,1]},{'tau~','tau'},[0 1]);
    tilde_f_u = NTset_block(tilde_f_u,{{'tau',1},{'tau~',1}},{[1,4],[3,2]},{'tau~','tau'},[-1; 0]);
    
    
 
    tilde_f_d = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[4,1],[2,2]},{'tau~','tau'},[-1; 0]);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[2,2],[1,4]},{'tau~','tau'},[0 1]);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[2,3],[1,1]},{'tau~','tau'},[1 0; 0 0; 0 1; 0 0]);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[1,1],[3,2]},{'tau~','tau'},[0 0 0 -1; 0 -1 0 0]);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[3,3],[5,1]},{'tau~','tau'},[0 1]);
    tilde_f_d = NTset_block(tilde_f_d,{{'tau',1},{'tau~',1}},{[1,5],[3,3]},{'tau~','tau'},[-1; 0]);
    
    
    %phase 
    ph = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[4,1],[4,1]},{'tau~','tau'}, [1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2,2],[2,2]},{'tau~','tau'}, -eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[2,3],[2,3]},{'tau~','tau'}, -eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,1],[1,1]},{'tau~','tau'}, eye(4));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3,2],[3,2]},{'tau~','tau'}, -eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[3,3],[3,3]},{'tau~','tau'}, -eye(2));
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[5,1],[5,1]},{'tau~','tau'}, [1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,5],[1,5]},{'tau~','tau'}, [1]);
    ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{[1,4],[1,4]},{'tau~','tau'}, [1]);
    
    %inf operator
    inf =  NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},2);
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[4,1],[4,1]},{'tau~','tau'}, [1]);
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[2,2],[2,2]},{'tau~','tau'}, -eye(2));
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[2,3],[2,3]},{'tau~','tau'}, -eye(2));
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[1,1],[1,1]},{'tau~','tau'}, eye(4));
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[3,2],[3,2]},{'tau~','tau'}, -eye(2));
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[3,3],[3,3]},{'tau~','tau'}, -eye(2));
    inf = NTset_block(inf,{{'tau',1},{'tau~',1}},{[5,1],[5,1]},{'tau~','tau'}, [1]);
    inf= NTset_block(inf,{{'tau',1},{'tau~',1}},{[1,5],[1,5]},{'tau~','tau'}, [1]);
    inf= NTset_block(inf,{{'tau',1},{'tau~',1}},{[1,4],[1,4]},{'tau~','tau'}, [1]);
    
    out = ABSITE_create(2,{ {[4,1],1},...
                            {[2,2],2},...
                            {[2,3],2},...
                            {[1,1],4},...
                            {[3,2],2},...
                            {[1,5],1},...
                            {[1,4],1},...
                            {[3,3],2},...
                            {[5,1],1}},'FERMIONS');

end



%% Construct the other operators and add them to the site.
% Spin-up occupation
n_u = NTdot(fdag_u, f_u, {'tau'},{'tau~'});
% Spin-down occupation
n_d = NTdot(fdag_d, f_d, {'tau'},{'tau~'});
%Total occupation
n = NTadd(n_u, n_d);

% Magnetic components.
Sz = NTadd(NTmult(NTdot(fdag_u, f_u, {'tau'},{'tau~'}), 0.5), NTmult(NTdot(fdag_d, f_d, {'tau'},{'tau~'}),-0.5 ));

Sp = NTdot(fdag_u, f_d, {'tau'},{'tau~'});
Sm = NTdot(fdag_d, f_u, {'tau'},{'tau~'});
%Hubbard coupling: fdag_u fdag_d f_d f_u
hub =  NTdot(n_u, n_d, {'tau'},{'tau~'});

% Superconducting components. 

Adag = NTdot(fdag_u, fdag_d, {'tau'},{'tau~'});
A = NTdot(f_d, f_u, {'tau'},{'tau~'});




% Spin-up occupation
tilde_n_u = NTdot(tilde_fdag_u, tilde_f_u, {'tau'},{'tau~'});
% Spin-down occupation
tilde_n_d = NTdot(tilde_fdag_d, tilde_f_d, {'tau'},{'tau~'});
%Total occupation
tilde_n = NTadd(tilde_n_u, tilde_n_d);

tilde_Sz = NTadd(NTmult(NTdot(tilde_fdag_u, tilde_f_u, {'tau'},{'tau~'}), 0.5), NTmult(NTdot(tilde_fdag_d, tilde_f_d, {'tau'},{'tau~'}),-0.5 ));
tilde_Sp = NTdot(tilde_fdag_u, tilde_f_d, {'tau'},{'tau~'});
tilde_Sm = NTdot(tilde_fdag_d, tilde_f_u, {'tau'},{'tau~'});
%Hubbard coupling: fdag_u fdag_d f_d f_u
tilde_hub = NTdot(tilde_n_u, tilde_n_d, {'tau'},{'tau~'});

Ddag_u =NTdot(fdag_u, tilde_fdag_u, {'tau'},{'tau~'});
Ddag_d =NTdot(fdag_d, tilde_fdag_d, {'tau'},{'tau~'});

Ddag = NTadd(Ddag_u, Ddag_d);

D_u = NTdot(f_u, tilde_f_u, {'tau'},{'tau~'});
D_d = NTdot(f_d, tilde_f_d, {'tau'},{'tau~'});

D = NTadd(D_u, D_d);

out.name = 'ABSITE_spinhalf_superfermion';
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
out = ABSITE_define_operator(out,'S+',Sp,1);
out = ABSITE_define_operator(out,'S-',Sm,1);


out = ABSITE_define_operator(out,'tilde_fdag_u',tilde_fdag_u,-1);
out = ABSITE_define_operator(out,'tilde_fdag_d',tilde_fdag_d,-1);
out = ABSITE_define_operator(out,'tilde_f_u',tilde_f_u,-1);
out = ABSITE_define_operator(out,'tilde_f_d',tilde_f_d,-1);
out = ABSITE_define_operator(out,'tilde_n_u',tilde_n_u,1);
out = ABSITE_define_operator(out,'tilde_n_d',tilde_n_d,1);
out = ABSITE_define_operator(out,'tilde_n',tilde_n,1);
out = ABSITE_define_operator(out,'tilde_hub',tilde_hub,1);

out = ABSITE_define_operator(out,'tilde_Sz',tilde_Sz,1);
out = ABSITE_define_operator(out,'tilde_S+',tilde_Sp,1);
out = ABSITE_define_operator(out,'tilde_S-',tilde_Sm,1);

out = ABSITE_define_operator(out,'Ddag_u',Ddag_u,1);
out = ABSITE_define_operator(out,'Ddag_d',Ddag_d,1);
out = ABSITE_define_operator(out,'D_u',D_u,1);
out = ABSITE_define_operator(out,'D_d',D_d,1);
out = ABSITE_define_operator(out,'D',D,1);
out = ABSITE_define_operator(out,'Ddag',Ddag,1);
out = ABSITE_define_operator(out,'Adag',Adag,1);
out = ABSITE_define_operator(out,'A',A,1);




end
