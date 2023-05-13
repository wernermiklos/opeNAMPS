function out = ABSITE_spinhalf( varargin )

%
%   Constructs a general site based on the U1_twoSz symmetry 
%
%
%   varargin can be a list with the names of the symmetries
%   Possible names to be used are  'U1_Q', 'Z1'
%   We can provide 0, 1 or 2 symmetries
%
%   Representation | 1   2  3   4  5   6  7  ...
%          2Sz     | 0  -1  1  -2  2  -3  3
%             
%
%  Possible states
%
%   States    2Sz    R_twoSz
%     |up>     1     3      
%     |dn>    -1     2     
%     
% . ↑ = char(8593)
% . ↓ = char(8595)
% 

if nargin > 1
     error('ABSITE_spinhalf() must have 1 argument');
end

% Read the names of the symmetries that we plan to use and their number
if nargin == 1
    symmetry_names = varargin{1};
    no_of_symmetries =numel(symmetry_names);
end

%if no arguments are provided that consider the U1_twoSz and the single symmetry.  
if nargin == 0
    symmetry_names={'U1_twoSz'};
    no_of_symmetries =numel(symmetry_names);
end

if  no_of_symmetries == 1 && strcmp(symmetry_names{1},'U1_twoSz')

Sz = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sz = NTset_block(Sz,{{'tau',1},{'tau~',1}},{[3],[3]},{'tau~','tau'},[0.5]);
Sz = NTset_block(Sz,{{'tau',1},{'tau~',1}},{[2],[2]},{'tau~','tau'},[-0.5]);

Sp = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sp = NTset_block(Sp,{{'tau',1},{'tau~',1}},{[2],[3]},{'tau~','tau'},[1]);

Sm = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sm = NTset_block(Sm,{{'tau',1},{'tau~',1}},{[3],[2]},{'tau~','tau'},[1]);

out = ABSITE_create(1,{{[2],1},{[3],1}},'BOSONS');

end

if   no_of_symmetries == 1 && strcmp(symmetry_names{1}, 'Z1')
    
    Sz = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    Sz = NTset_block(Sz,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [ 0.5 0; 0 -0.5]);
    
    Sp = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    Sp = NTset_block(Sp,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 1; 0 0]);
    
    Sm = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
    Sm = NTset_block(Sm,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'}, [0 0; 1 0]);
        
    out = ABSITE_create(1,{{[1],2}},'BOSONS');
    
end

% Sx and Sy components of the spin

Sx = NTadd(NTmult(Sp, 0.5 ), NTmult(Sm, 0.5 ));
Sy = NTadd(NTmult(Sp, -0.5i ), NTmult(Sm,0.5i));

out.name = 'ABSITE_spinhalf';

out = ABSITE_define_operator(out,'Sz',Sz,1);
out = ABSITE_define_operator(out,'S+',Sp,1);
out = ABSITE_define_operator(out,'S-',Sm,1);
out = ABSITE_define_operator(out,'Sx',Sx,1);
out = ABSITE_define_operator(out,'Sy',Sy,1);

end