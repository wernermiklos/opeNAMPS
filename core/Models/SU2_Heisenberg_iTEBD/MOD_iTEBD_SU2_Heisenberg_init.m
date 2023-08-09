function TEBD_env = MOD_iTEBD_SU2_Heisenberg_init(SU2sym,unitcell_length,spinsize, J,dt)
%MOD_TEBD_SU2XU1_HUBBARD_INIT Summary of this function goes here
%   Detailed explanation goes here
    NO_OF_SYMMETRIES = 1;
    
    % TEBD environment initialization:
    EmptyMPS = LSMPS_create(unitcell_length,NO_OF_SYMMETRIES);
    spin_site = SITE_generate_SU2_spin_site(SU2sym,spinsize);
    Sites = cellfun(@(x) spin_site,cell(1,unitcell_length),'UniformOutput',false);  % cell array of spinhalf sites
    TEBD_env = TEBD_init(EmptyMPS,{SU2sym},Sites,'INFINITE',true);

    % Building up the Heisenberg Hamiltonian
    % the spin tensor operators components are {-sqrt(0.5)*Sp, Sz, sqrt(0.5)*Sm}, and
    % the Heisenberg coupling is (0.5*SpSm + SzSz + 0.5*SmSp)
    heis_coup = COUP_generate_TwoSite_full_OpProd({SU2sym},spin_site,spin_site,'S','S',...
                                                 {{[1,3],-1},{[2,2],1},{[3,1],-1}}); 

    for bond_pos = 1:unitcell_length
        TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env,bond_pos,{{heis_coup,J}});
    end

    % Determination of two site evolvers (both full and reduced);
    TEBD_env = TEBD_generate_TwoSite_full_evolver_list(TEBD_env,{SU2sym},dt);
    TEBD_env = TEBD_generate_reduced_TwoSite_ops(TEBD_env);
end

