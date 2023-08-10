% Non-Abelian (SU(2)) iTEBD script for the Heisenberg chain.
 S = 2.0;    % size of the single spin

 locrep_max = (2*S+1)*2;     % The maximal local representation (for larger spins larger values are needed...)
 bondrep_max = 3*locrep_max;   % The maximal bond representation

% Generating symmetries
 SU2sym = generate_SU2_Symmetry(1:bondrep_max,1:locrep_max,1:bondrep_max);

% Model parameters 
 L = 20;     % chain length
 J = 1.0;    % the coupling strength
 dt=0.01;    % time step length
 Nstep = 500;  % number of time steps
 M_mult = 512; % max (multiplet) bond dimension
 
% Charge conventions:  (Q_empty = -1, Q_fermion = 1 is suitable at half
% filling)
Q_empty = -1;
Q_fermion = 1;

rep_site = 2*S+1;
rep_triv = 1;


% Creation of TEBD environment.
TEBD_env = MOD_TEBD_SU2_Heisenberg_init(SU2sym,...
                                         L, ...
                                         S, ...     %spin size
                                         J, ...     %coupling strength
                                         dt);       %time step


% Setting up the initial state (singlet bond at every second bond)
for pos = 1:L
    if mod(pos,2)
        TEBD_env.StateMPS.matrices{pos} = NTset_block(TEBD_env.StateMPS.matrices{pos},...
                                                   {{'t_left',1},{'tau',1},{'t_right',1}},...
                                                   {rep_triv,rep_site,rep_site},...
                                                   {'t_left','tau','alpha','t_right'},...
                                                   [1]);
        TEBD_env.StateMPS.schmidt_list{pos} = NTset_block(TEBD_env.StateMPS.schmidt_list{pos}, ...
                                                        {{'t_left',1}}, ...
                                                        {rep_site},...
                                                        {'t_left','t_right'}, ...
                                                        [1]);
    else
        TEBD_env.StateMPS.matrices{pos} = NTset_block(TEBD_env.StateMPS.matrices{pos},...
                                                   {{'t_left',1},{'tau',1},{'t_right',1}},...
                                                   {rep_site,rep_site,rep_triv},...
                                                   {'t_left','tau','alpha','t_right'},...
                                                   [1]);
        TEBD_env.StateMPS.schmidt_list{pos} = NTset_block(TEBD_env.StateMPS.schmidt_list{pos}, ...
                                                        {{'t_left',1}}, ...
                                                        {rep_triv},...
                                                        {'t_left','t_right'}, ...
                                                [1]);
    end
end

% Simulate dynamics
bond_energy = zeros(0,L-1);   % expectation value of bond energy will be stored here
bond_entropy = zeros(0,L-1);   % expectation value of bond energy will be stored here 
t = (0:dt:(Nstep*dt))';

% Measurement of average charge density on both sublattices
for pos = 1:(L-1)
    bond_energy(1,pos) = real(LSMPS_expval_scalar_bond_operator(TEBD_env.StateMPS,...
                                                                pos,...
                                                                TEBD_env.TwoSiteHlist_red{pos}));
end

for i = 1:Nstep
    % 
    [TEBD_env,info] = TEBD_evolve_Trotter1_Ured(TEBD_env,M_mult);  % here precalculated reduced evolvers are used 
    %[TEBD_env,info] = TEBD_evolve_Trotter1_Ufull(TEBD_env,M_mult);   % here full evolvers are used, and Clebsch-factor is used on the fly.


    % Measurement of average charge density on both sublattices 
    for pos = 1:(L-1)
        bond_energy(i+1,pos) = real(LSMPS_expval_scalar_bond_operator(TEBD_env.StateMPS,...
                                                                      pos,...
                                                                      TEBD_env.TwoSiteHlist_red{pos}));
        bond_entropy(i+1,pos) = LSMPS_measure_bond_entropy(TEBD_env.StateMPS,...
                                                                {SU2sym}, ...
                                                                pos);
    end
    disp(['time= ', num2str(i*dt)]);
    disp(['Bond dimensions: ', num2str(info.M)]);
    disp(['Bond energies: ', num2str(bond_energy(i+1,:))]);
    disp(['Bond entropies: ', num2str(bond_entropy(i+1,:))]);
end

