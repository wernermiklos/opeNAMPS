function TEBD_env = MOD_iTEBD_U1_NN_spinless_fermion_init(U1sym,unitcell_length,V,dt,varargin)
%MOD_TEBD_SU2XU1_HUBBARD_INIT Summary of this function goes here
%   Detailed explanation goes here
    NO_OF_SYMMETRIES = 1;

    parameters = struct();
       % Default parameters
    parameters.Q_EMPTY = -1;
    parameters.Q_OCCUPIED = +1;
       % Update according to varargin
    parameters = parameter_updater(parameters,varargin);
   
    
    % TEBD environment initialization:
    EmptyMPS = LSMPS_create(unitcell_length,NO_OF_SYMMETRIES);
    spinless_site = SITE_generate_U1_spinless_fermion(U1sym,...
                                                      parameters.Q_EMPTY, ...
                                                      parameters.Q_OCCUPIED);
    Sites = cellfun(@(x) spinless_site,cell(1,unitcell_length),'UniformOutput',false);  % cell array of spinhalf sites
    TEBD_env = TEBD_init(EmptyMPS,{U1sym},Sites,'INFINITE',true,'FERMIONS',true);

    % Building up the Hamiltonian
    hopp_lr = COUP_generate_TwoSite_full_OpProd({U1sym},spinless_site,spinless_site,'fdag','f',...
                                                 {{[1,1],-1}});
    hopp_rl = COUP_generate_TwoSite_full_OpProd({U1sym},spinless_site,spinless_site,'f','fdag',...
                                                 {{[1,1],1}});
    interact = COUP_generate_TwoSite_full_OpProd({U1sym},spinless_site,spinless_site,'n','n',...
                                                 {{[1,1],V}});
    for bond_pos = 1:unitcell_length
        TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env,bond_pos,{{hopp_lr,1},{hopp_rl,1},{interact,1}});
    end

    % Determination of two site evolvers (both full and reduced);
    TEBD_env = TEBD_generate_TwoSite_full_evolver_list(TEBD_env,{U1sym},dt);
    TEBD_env = TEBD_generate_reduced_TwoSite_ops(TEBD_env);
end

