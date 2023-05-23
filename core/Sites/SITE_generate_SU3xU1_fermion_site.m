function SU3site = SITE_generate_SU3xU1_fermion_site(SU3sym, U1sym)
%GENERATE_SU3_SITE Generates a three-color SU(3) symmetric fermionic site
%   Detailed explanation goes here
SU3site = SITE_create({SU3sym,U1sym},{{[1,2],1},{[2,1],1},{[3,3],1},{[1,5],1}});
fdag_list = cell(1,3);
f_list = cell(1,3);

fdag_list{1} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
fdag_list{1} = NTset_block(fdag_list{1}, {{'tau~',1},{'tau',1}},{[2,1],[1,2]},{'mu~','mu','tau~','tau'},reshape([1; 0; 0],[3,1,1,1]));
fdag_list{1} = NTset_block(fdag_list{1}, {{'tau~',1},{'tau',1}},{[3,3],[2,1]},{'mu~','mu','tau~','tau'},reshape([0, 1, 0; ...
                                                                                               0, 0, 1; ...
                                                                                               0, 0, 0], [3,3,1,1]));
fdag_list{1}= NTset_block(fdag_list{1}, {{'tau~',1},{'tau',1}},{[1,5],[3,3]},{'mu~','mu','tau~','tau'},reshape([0, 0, 1],[1,3,1,1]));


f_list{1} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
f_list{1} = NTset_block(f_list{1}, {{'tau~',1},{'tau',1}},{[1,2],[2,1]},{'mu~','mu','tau~','tau'},reshape([0, 0, 1],[1,3,1,1]));
f_list{1} = NTset_block(f_list{1}, {{'tau~',1},{'tau',1}},{[2,1],[3,3]},{'mu~','mu','tau~','tau'},reshape([0, -1, 0; ...
                                                                                            0, 0, -1; ...
                                                                                            0, 0, 0], [3,3,1,1]));
f_list{1} = NTset_block(f_list{1}, {{'tau~',1},{'tau',1}},{[3,3],[1,5]},{'mu~','mu','tau~','tau'},reshape([1; 0; 0],[3,1,1,1]));


fdag_list{2} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
fdag_list{2} = NTset_block(fdag_list{2}, {{'tau~',1},{'tau',1}},{[2,1],[1,2]},{'mu~','mu','tau~','tau'},reshape([0; 1; 0],[3,1,1,1]));
fdag_list{2} = NTset_block(fdag_list{2}, {{'tau~',1},{'tau',1}},{[3,3],[2,1]},{'mu~','mu','tau~','tau'},reshape([-1, 0, 0; ...
                                                                                               0, 0, 0; ...
                                                                                               0, 0, 1], [3,3,1,1]));
fdag_list{2} = NTset_block(fdag_list{2}, {{'tau~',1},{'tau',1}},{[1,5],[3,3]},{'mu~','mu','tau~','tau'},reshape([0, -1, 0],[1,3,1,1]));


f_list{2} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
f_list{2} = NTset_block(f_list{2}, {{'tau~',1},{'tau',1}},{[1,2],[2,1]},{'mu~','mu','tau~','tau'},reshape([0, -1, 0],[1,3,1,1]));
f_list{2} = NTset_block(f_list{2}, {{'tau~',1},{'tau',1}},{[2,1],[3,3]},{'mu~','mu','tau~','tau'},reshape([1, 0, 0; ...
                                                                                            0, 0, 0; ...
                                                                                            0, 0, -1], [3,3,1,1]));
f_list{2} = NTset_block(f_list{2}, {{'tau~',1},{'tau',1}},{[3,3],[1,5]},{'mu~','mu','tau~','tau'},reshape([0; 1; 0],[3,1,1,1]));


fdag_list{3} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
fdag_list{3} = NTset_block(fdag_list{3}, {{'tau~',1},{'tau',1}},{[2,1],[1,2]},{'mu~','mu','tau~','tau'},reshape([0; 0; 1],[3,1,1,1]));
fdag_list{3} = NTset_block(fdag_list{3}, {{'tau~',1},{'tau',1}},{[3,3],[2,1]},{'mu~','mu','tau~','tau'},reshape([0, 0, 0; ...
                                                                                               -1, 0, 0; ...
                                                                                               0, -1, 0], [3,3,1,1]));
fdag_list{3} = NTset_block(fdag_list{3}, {{'tau~',1},{'tau',1}},{[1,5],[3,3]},{'mu~','mu','tau~','tau'},reshape([1, 0, 0],[1,3,1,1]));

f_list{3} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
f_list{3} = NTset_block(f_list{3}, {{'tau~',1},{'tau',1}},{[1,2],[2,1]},{'mu~','mu','tau~','tau'},reshape([1, 0, 0],[1,3,1,1]));
f_list{3} = NTset_block(f_list{3}, {{'tau~',1},{'tau',1}},{[2,1],[3,3]},{'mu~','mu','tau~','tau'},reshape([0, 0, 0; ...
                                                                                            1, 0, 0; ...
                                                                                            0, 1, 0], [3,3,1,1]));
f_list{3} = NTset_block(f_list{3}, {{'tau~',1},{'tau',1}},{[3,3],[1,5]},{'mu~','mu','tau~','tau'},reshape([0; 0; 1],[3,1,1,1]));


ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
ph = NTset_block(ph, {{'tau~',1},{'tau',1}},{[1,2],[1,2]},{'mu~','mu','tau~','tau'},reshape([1],[1,1,1,1]));
ph = NTset_block(ph, {{'tau~',1},{'tau',1}},{[2,1],[2,1]},{'mu~','mu','tau~','tau'},reshape(-eye(3),[3,3,1,1]));
ph = NTset_block(ph, {{'tau~',1},{'tau',1}},{[3,3],[3,3]},{'mu~','mu','tau~','tau'},reshape(eye(3),[3,3,1,1]));
ph = NTset_block(ph, {{'tau~',1},{'tau',1}},{[1,5],[1,5]},{'mu~','mu','tau~','tau'},reshape([-1],[1,1,1,1]));


n = NTadd(NTdot(fdag_list{1},f_list{3}, {'tau', 'mu'}, {'tau~', 'mu~'}), ...
          NTadd(NTneg(NTdot(fdag_list{2},f_list{2}, {'tau', 'mu'}, {'tau~', 'mu~'})), ...
                NTdot(fdag_list{3},f_list{1}, {'tau', 'mu'}, {'tau~', 'mu~'})));
hubbard_cup = NTmult(NTsubtr(NTdot(n, n, {'tau', 'mu'}, {'tau~', 'mu~'}),n), 0.5);

SU3site = SITE_define_tensor_operator(SU3site, 'fdag',{SU3sym,U1sym},[2,3],fdag_list,-1);
SU3site = SITE_define_tensor_operator(SU3site, 'f',{SU3sym,U1sym},[3,2],f_list,-1);
SU3site = SITE_define_tensor_operator(SU3site, 'ph',{SU3sym,U1sym},[1,1],{ph},1);
SU3site = SITE_define_tensor_operator(SU3site, 'n',{SU3sym,U1sym},[1,1],{n},1);
SU3site = SITE_define_tensor_operator(SU3site, 'hub',{SU3sym,U1sym},[1,1],{hubbard_cup},1);


end

