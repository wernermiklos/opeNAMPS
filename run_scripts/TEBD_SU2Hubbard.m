% Non-Abelian (SU(2) x U(1)) iTEBD script for the Hubbard chain.
 
% Generating symmetries
 SU2sym = generate_SU2_Symmetry(1:31,1:11,1:41);
 U1sym = generate_U1_Symmetry(1:31,1:11,1:41);

% Model parameters 
 unitcell_length = 2;
 U=1.0;
 dt=0.01;
 Nstep = 500;
 M_mult = 1024;
 
% Charge conventions:  (Q_empty = -1, Q_fermion = 1 is suitable at half
% filling)
Q_empty = -1;
Q_fermion = 1;

rep_empty = [SU2sym.Qnum_to_Gamma(0), U1sym.Qnum_to_Gamma(Q_empty)];
rep_single = [SU2sym.Qnum_to_Gamma(0.5) ,U1sym.Qnum_to_Gamma(Q_empty + Q_fermion)];
rep_double = [SU2sym.Qnum_to_Gamma(0), U1sym.Qnum_to_Gamma(Q_empty + 2*Q_fermion)];
rep_triv = [1,1];


% Creation of TEBD environment.
TEBD_env = MOD_iTEBD_SU2xU1_spinhalf_Hubbard_init(SU2sym,U1sym,...
                                                  unitcell_length, ...
                                                  U, ...
                                                  dt, ...
                                                  'CHARGE_OF_EMPTY_SITE',Q_empty, ...
                                                  'CHARGE_PER_FERMION',Q_fermion);


% Setting up the initial state ( (updown) singlets at every second site)
TEBD_env.StateMPS.matrices{1} = NTset_block(TEBD_env.StateMPS.matrices{1},...
                                           {{'t_left',1},{'tau',1},{'t_right',1}},...
                                           {rep_triv,rep_double,rep_double},...
                                           {'t_left','tau','alpha','t_right'},...
                                           [1]);
TEBD_env.StateMPS.schmidt_list{1} = NTset_block(TEBD_env.StateMPS.schmidt_list{1}, ...
                                                {{'t_left',1}}, ...
                                                {rep_double},...
                                                {'t_left','t_right'}, ...
                                                [1]);
TEBD_env.StateMPS.matrices{2} = NTset_block(TEBD_env.StateMPS.matrices{2},...
                                           {{'t_left',1},{'tau',1},{'t_right',1}},...
                                           {rep_double,rep_empty,rep_triv},...
                                           {'t_left','tau','alpha','t_right'},...
                                           [1]);
TEBD_env.StateMPS.schmidt_list{2} = NTset_block(TEBD_env.StateMPS.schmidt_list{2}, ...
                                                {{'t_left',1}}, ...
                                                {rep_triv},...
                                                {'t_left','t_right'}, ...
                                                [1]);

% Simulate dynamics
n_avg = zeros(0,2);
t = (0:dt:(Nstep*dt))';
n_avg(1,1) = LSMPS_expval_local_scalar_operator(TEBD_env.StateMPS,...
                                                1,...
                                                TEBD_env.Sites{1}.operators('n'));
n_avg(1,2) = LSMPS_expval_local_scalar_operator(TEBD_env.StateMPS,...
                                                2,...
                                                TEBD_env.Sites{2}.operators('n'));
for i = 1:Nstep
    [TEBD_env,info] = TEBD_evolve_Trotter1_Ured(TEBD_env,M_mult);
    n_avg(end+1,1) = LSMPS_expval_local_scalar_operator(TEBD_env.StateMPS,...
                                                       1,...
                                                       TEBD_env.Sites{1}.operators('n'));
    n_avg(end,2) = LSMPS_expval_local_scalar_operator(TEBD_env.StateMPS,...
                                                     2,...
                                                     TEBD_env.Sites{2}.operators('n'));
    disp([i*dt,n_avg(end,1),n_avg(end,2),info.M(1),info.M(2)]);
end

% a comment
