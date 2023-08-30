% Non-Abelian (SU(2)) iTEBD script for the Heisenberg chain.
 S = 2.0;    % size of the single spin

 locrep_max = (2*S+1)*2;     % The maximal local representation (for larger spins larger values are needed...)
 bondrep_max = 3*locrep_max;   % The maximal bond representation

% Generating symmetries
 SU2sym = generate_SU2_Symmetry(1:bondrep_max,1:locrep_max,1:bondrep_max);

% Model parameters 
 unitcell_length = 2;
 J = 1.0;    %the coupling strength
 dt=0.01;
 Nstep = 500;
 M_mult = 512;
 
% Charge conventions:  (Q_empty = -1, Q_fermion = 1 is suitable at half
% filling)
%Q_empty = -1;
%Q_fermion = 1;

rep_site = 2*S+1;
rep_triv = 1;


% Creation of TEBD environment.
TEBD_env = MOD_iTEBD_SU2_Heisenberg_init(SU2sym,...
                                         unitcell_length, ...
                                         S, ...     %spin size
                                         J, ...     %coupling strength
                                         dt);       %time step


% Setting up the initial state (singlet bond at every second bond)
TEBD_env.StateMPS.matrices{1} = NTset_block(TEBD_env.StateMPS.matrices{1},...
                                           {{'t_left',1},{'tau',1},{'t_right',1}},...
                                           {rep_triv,rep_site,rep_site},...
                                           {'t_left','tau','alpha','t_right'},...
                                           [1]);
TEBD_env.StateMPS.schmidt_list{1} = NTset_block(TEBD_env.StateMPS.schmidt_list{1}, ...
                                                {{'t_left',1}}, ...
                                                {rep_site},...
                                                {'t_left','t_right'}, ...
                                                [1]);
TEBD_env.StateMPS.matrices{2} = NTset_block(TEBD_env.StateMPS.matrices{2},...
                                           {{'t_left',1},{'tau',1},{'t_right',1}},...
                                           {rep_site,rep_site,rep_triv},...
                                           {'t_left','tau','alpha','t_right'},...
                                           [1]);
TEBD_env.StateMPS.schmidt_list{2} = NTset_block(TEBD_env.StateMPS.schmidt_list{2}, ...
                                                {{'t_left',1}}, ...
                                                {rep_triv},...
                                                {'t_left','t_right'}, ...
                                                [1]);

% Simulate dynamics
bond_energy = zeros(0,2);   % expectation value of bond energy will be stored here 
t = (0:dt:(Nstep*dt))';

% Measurement of average charge density on both sublattices
for pos = 1:2
    bond_energy(1,pos) = real(LSMPS_expval_scalar_bond_operator(TEBD_env.StateMPS,...
                                                                pos,...
                                                                TEBD_env.TwoSiteHlist_red{pos}));
end

for i = 1:Nstep
    % 
    [TEBD_env,info] = TEBD_evolve_Trotter1_Ured(TEBD_env,M_mult);  % here precalculated reduced evolvers are used 
    %[TEBD_env,info] = TEBD_evolve_Trotter1_Ufull(TEBD_env,M_mult);   % here full evolvers are used, and Clebsch-factor is used on the fly.


    % Measurement of average charge density on both sublattices 
    for pos = 1:2
        bond_energy(i+1,pos) = real(LSMPS_expval_scalar_bond_operator(TEBD_env.StateMPS,...
                                                                      pos,...
                                                                      TEBD_env.TwoSiteHlist_red{pos}));
    end
    disp([i*dt,bond_energy(end,1),bond_energy(end,2),info.M(1),info.M(2)]);
end

% a comment

