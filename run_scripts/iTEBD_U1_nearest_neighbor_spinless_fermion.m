dt = 0.01;  % time step
V = 2;      % nearest neighbor interaction strength
unitcell_length=2;        % unit cell length
NStep = 100; 
Mmax = 64; 


NO_OF_SYMMETRIES = 1;     % Do not change, unless you rewrite the script for more symmetries. 

U1Sym = generate_U1_Symmetry(1:41,1:21,1:61);    % in the non-Abelian code we have to set the 
                                                 % possible representation indices.

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


TEBD_env = MOD_iTEBD_U1_NN_spinless_fermion_init(U1sym,unitcell_length,V,dt,'Q_EMPTY',Q_EMPTY,'Q_OCCUPIED',Q_OCCUPIED); 

for pos=1:unitcell_length
    switch mod(pos,2)
        case 1
            TEBD_env.StateMPS.matrices{pos} = NTset_block(TEBD_env.StateMPS.matrices{pos}, ...
                                                 {{'t_left',1}, ...
                                                  {'tau',1}, ...
                                                  {'t_right',1}}, ...
                                                  {REP_TRIV,REP_OCCUPIED,REP_OCCUPIED}, ...
                                                  {'t_left','tau','alpha','t_right'}, ...
                                                  1);
            TEBD_env.StateMPS.schmidt_list{pos} = NTset_block(TEBD_env.StateMPS.schmidt_list{pos},...
                                                     {{'t_left',1}}, ...
                                                     {REP_OCCUPIED}, ...
                                                     {'t_left','t_right'}, ...
                                                     1);
        case 0
            TEBD_env.StateMPS.matrices{pos} = NTset_block(TEBD_env.StateMPS.matrices{pos}, ...
                                                 {{'t_left',1}, ...
                                                  {'tau',1}, ...
                                                  {'t_right',1}}, ...
                                                  {REP_OCCUPIED,REP_EMPTY,REP_TRIV}, ...
                                                  {'t_left','tau','alpha','t_right'}, ...
                                                  1);
            TEBD_env.StateMPS.schmidt_list{pos} = NTset_block(TEBD_env.StateMPS.schmidt_list{pos},...
                                                     {{'t_left',1}}, ...
                                                     {REP_TRIV}, ...
                                                     {'t_left','t_right'}, ...
                                                     1);
    end
end



n_avg = internal_measure_singlesite_charges(TEBD_env);
t = 0;
disp(['time = ', num2str(t)]);
disp(['bond_dims = ', num2str(ones(1,unitcell_length))]);
disp(['log10(trunc_errs) = []']);
disp(['<n_i> = ', num2str(n_avg(end,:))]);



for i = 1:NStep 
    tic;
    [TEBD_env, info] = TEBD_evolve_Trotter1_Ured(TEBD_env, Mmax);
    %measurement of n(x=1)
    n_avg = [n_avg; internal_measure_singlesite_charges(TEBD_env)];
    disp(['time = ', num2str(i*dt)]);
    disp(['bond_dims = ', num2str(info.M)]);
    disp(['log10(truncerr) = ', num2str(real(cellfun(@log10, info.truncerr)))]);
    disp(['<n_i> = ', num2str(n_avg(end,:))]);
    toc;
end

function n_avg_snapshot = internal_measure_singlesite_charges(TEBD_env)
   n_avg_snapshot = zeros(1,TEBD_env.chain_length);
    for pos = 1:TEBD_env.chain_length
       tmp_mx = NTdot(TEBD_env.StateMPS.matrices{pos}, TEBD_env.StateMPS.schmidt_list{pos}, {'t_right'},{'t_left'});
       tmp = NTdot(tmp_mx,TEBD_env.Sites{pos}.operators('n'),{'tau'},{'tau'});
       n_avg_snapshot(pos) = NTdot(tmp,NTconj(tmp_mx),{'t_left','tau~','alpha','t_right'},{'t_left','tau','alpha','t_right'});
    end
end