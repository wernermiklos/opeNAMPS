dt = 0.01;  % time step
V = 2;      % nearest neighbor interaction strength
L=2;        % unit cell length, if INFINITE is true, chain length otherwise.
NStep = 100; 
Mmax = 64; 

NO_OF_SYMMETRIES = 1;     % Do not change, unless you rewrite the script for more symmetries.
INFINITE = true;   

U1Sym = generate_U1_Symmetry([],[],[]);    % in the Abelian code it is enough to generate the symm without the Clebsch tensor. 

% In the infinite case filling is essential: the total 'charge' of the unit
% cell must be zero. I.e. we need to specify the charge of the empty and
% occupied sites so, that the total charge will be zero.
% E.g. for half filling Q_empty = -1 and Q_occupied = +1 will work.
% For third filling Q_empty = -1 and Q_occupied = +2 is the good choice.

%half filling:
Q_EMPTY = -1;
Q_OCCUPIED = +1;
REP_EMPTY = U1Sym.Qnum_to_Gamma(Q_EMPTY);
REP_OCCUPIED = U1Sym.Qnum_to_Gamma(Q_OCCUPIED);
REP_TRIV = 1;
% Initial state: CDW state (every second site is occupied: 10101010);

StartMPS = ABMPS_create(L,L,NO_OF_SYMMETRIES);    % cut is at the end of the chain! (left canonical MPS)
for pos=1:L
    switch mod(pos,2)
        case 1
            A = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},NO_OF_SYMMETRIES); 
            A = NTset_block(A,{{'t_in',1},{'tau',1},{'t_out',1}},{REP_TRIV,REP_OCCUPIED,REP_OCCUPIED},{'t_in','tau','t_out'},[1]); 
            StartMPS = ABMPS_set_matrix(StartMPS, pos, A); 
            schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},NO_OF_SYMMETRIES); 
            schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{REP_OCCUPIED,U1Sym.conj_reps(REP_OCCUPIED)},{'t_left','t_right'},[1]);
            StartMPS.schmidt_list{pos} = schmidt; 
        case 0
            A = NAtensor({'t_in','tau','t_out'},{'i','i','o'},{[1],[2],[3]},NO_OF_SYMMETRIES); 
            A = NTset_block(A,{{'t_in',1},{'tau',1},{'t_out',1}},{REP_OCCUPIED, REP_EMPTY, REP_TRIV},{'t_in','tau','t_out'},[1]); 
            StartMPS = ABMPS_set_matrix(StartMPS, pos, A); 
            schmidt = NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},NO_OF_SYMMETRIES); 
            schmidt = NTset_block(schmidt,{{'t_left',1},{'t_right',1}},{REP_TRIV,U1Sym.conj_reps(REP_TRIV)},{'t_left','t_right'},[1]);   % remark: conj rep of triv is triv.
            StartMPS.schmidt_list{pos} = schmidt; 
    end
end


StartMPS = ABMPS_set_schmidt(StartMPS,StartMPS.schmidt_list{end}); 

ABTEBD_env = MOD_ABTEBD_U1_nearest_neighbor_spinless_fermion_init(U1sym,StartMPS,V,dt,'INFINITE',INFINITE,'Q_EMPTY',Q_EMPTY,'Q_OCCUPIED',Q_OCCUPIED); 

n_avg = internal_measure_singlesite_charges(ABTEBD_env);
t = 0;
disp(['time = ', num2str(t)]);
disp(['bond_dims = ', num2str(ones(1,L))]);
disp(['log10(trunc_errs) = []']);
disp(['<n_i> = ', num2str(n_avg(end,:))]);





for i = 1:NStep 
    tic;
    [ABTEBD_env, info] = ABTEBD_evolve_Trotter1(ABTEBD_env, Mmax, 1);
    %measurement of n(x=1)
    n_avg = [n_avg; internal_measure_singlesite_charges(ABTEBD_env)];
    disp(['time = ', num2str(i*dt)]);
    disp(['bond_dims = ', num2str(info.bond_dims)]);
    disp(['log10(trunc_errs) = ', num2str(real(log10(info.trunc_errs)))]);
    disp(['<n_i> = ', num2str(n_avg(end,:))]);
    toc;
end

function n_avg_snapshot = internal_measure_singlesite_charges(ABTEBD_env)
   n_avg_snapshot = zeros(1,ABTEBD_env.chain_length);
    for pos = 1:ABTEBD_env.chain_length
       tmp_mx = NTdot(ABTEBD_env.StateMPS.left_matrices{pos}, ABTEBD_env.StateMPS.schmidt_list{pos}, {'t_out'},{'t_left'},{{},{{'t_right','t_out'}}});
       tmp = NTdot(tmp_mx,ABTEBD_env.Sites{pos}.operators('n'),{'tau'},{'tau'});
       n_avg_snapshot(pos) = NTdot(tmp,NTconj(tmp_mx),{'t_in','tau~','t_out'},{'t_in','tau','t_out'});
    end
end