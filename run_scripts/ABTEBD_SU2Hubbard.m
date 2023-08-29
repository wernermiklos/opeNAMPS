dt = 0.01; 
U = 3;
L=16;
NStep = 100; 
Mmax = 64; 

U1Sym = generate_U1_Symmetry(1:31,1:21,1:51); 

%first U1: spin-Z   (up: rep = 3, down: rep = 2)
%second U1:  charge
%empty site:      charge -1  (rep = 2)
%single occ site: charge 0   (rep = 1)
%double occ site: charge 1   (rep = 3)
%initial state:  CDW with double occ sites at every second position

%U1sym = generate_U1_Symmetry(1:31,1:21,1:51);

StartMPS = ABMPS_create(L,L,2); 
for pos=1:L/2
    sitepos=2*pos-1;
    A = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},2); 
    A = NTset_block(A,{{'t_in',1},{'tau',1},{'t_out',1}},{[1,1],[1,3],[1,3]},{'t_in','tau','t_out'},[1]); 
    StartMPS = ABMPS_set_matrix(StartMPS,sitepos,A); 
    schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},2); 
    schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{[1,3],SUPER_conj_rep({U1Sym,U1Sym},[1,3])},{'t_left','t_right'},[1]);
    StartMPS.schmidt_list{sitepos} = schmidt; 
    
    sitepos=2*pos;
    A = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},2); 
    A = NTset_block(A,{{'t_in',1},{'tau',1},{'t_out',1}},{[1,3],[1,2],[1,1]},{'t_in','tau','t_out'},[1]); 
    StartMPS = ABMPS_set_matrix(StartMPS,sitepos,A); 
    schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},2); 
    schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{[1,1],SUPER_conj_rep({U1Sym,U1Sym},[1,1])},{'t_left','t_right'},[1]);
    StartMPS.schmidt_list{sitepos} = schmidt; 
end


StartMPS = ABMPS_set_schmidt(StartMPS,schmidt); 

ABTEBD_env = MOD_ABTEBD_spinhalf_hubbard_init(U1Sym,StartMPS, U, dt); 

n_avg = internal_measure_singlesite_charges(ABTEBD_env);
t = 0;
disp(['time = ', num2str(t)]);
disp(['bond_dims = ', num2str(ones(1,L))]);
disp(['log10(trunc_errs) = []']);
disp(['<n_i> = ', num2str(n_avg(end,:))]);
% tmp_mx = NTdot(ABTEBD_env.StateMPS.left_matrices{1}, ABTEBD_env.StateMPS.schmidt_list{1},{'t_out'},{'t_left'},{{},{{'t_right','t_out'}}});
% tmp = NTdot(tmp_mx,ABTEBD_env.Sites{1}.operators('n'),{'tau'},{'tau'});
% NTdot(tmp,NTconj(tmp_mx),{'t_in','tau~','t_out'},{'t_in','tau','t_out'})


for i = 1:NStep 
    tic;
    [ABTEBD_env, info] = ABTEBD_evolve_Trotter1(ABTEBD_env, Mmax, 1);
    %measurement of n(x=1)
%     tmp_mx = NTdot(ABTEBD_env.StateMPS.left_matrices{1}, ABTEBD_env.StateMPS.schmidt_list{1},{'t_out'},{'t_left'},{{},{{'t_right','t_out'}}});
%     tmp = NTdot(tmp_mx,ABTEBD_env.Sites{1}.operators('n'),{'tau'},{'tau'});
%     navg(end+1) = NTdot(tmp,NTconj(tmp_mx),{'t_in','tau~','t_out'},{'t_in','tau','t_out'});
%     tmp_mx = NTdot(ABTEBD_env.StateMPS.left_matrices{2}, ABTEBD_env.StateMPS.schmidt_list{2},{'t_out'},{'t_left'},{{},{{'t_right','t_out'}}});
%     tmp = NTdot(tmp_mx,ABTEBD_env.Sites{2}.operators('n'),{'tau'},{'tau'});
%     navg2(end+1) = NTdot(tmp,NTconj(tmp_mx),{'t_in','tau~','t_out'},{'t_in','tau','t_out'});
%     disp(real([i*dt, navg(end), navg2(end), info.bond_dims(1), info.bond_dims(2), real(log10(info.trunc_errs(1))), real(log10(info.trunc_errs(2)))]));
    n_avg = [n_avg; internal_measure_singlesite_charges(ABTEBD_env)];
    disp(['time = ', num2str(i*dt)]);
    disp(['bond_dims = ', num2str(info.bond_dims)]);
    disp(['log10(trunc_errs) = ', num2str(real(log10(info.trunc_errs)))]);
    disp(['<n_i> = ', num2str(n_avg(end,:))]);
    toc;
end
%plot(dt*(1:NStep),navg)

function n_avg_snapshot = internal_measure_singlesite_charges(ABTEBD_env)
   n_avg_snapshot = zeros(1,ABTEBD_env.chain_length);
    for pos = 1:ABTEBD_env.chain_length
       tmp_mx = NTdot(ABTEBD_env.StateMPS.left_matrices{pos}, ABTEBD_env.StateMPS.schmidt_list{pos}, {'t_out'},{'t_left'},{{},{{'t_right','t_out'}}});
       tmp = NTdot(tmp_mx,ABTEBD_env.Sites{pos}.operators('n'),{'tau'},{'tau'});
       n_avg_snapshot(pos) = NTdot(tmp,NTconj(tmp_mx),{'t_in','tau~','t_out'},{'t_in','tau','t_out'});
    end
end