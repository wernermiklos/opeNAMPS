function TEBD_env = MOD_TEBD_SU2xU1_spinhalf_Hubbard_open_init(SU2sym,U1sym,L,U,G,dt)
%MOD_TEBD_SU2XU1_HUBBARD_INIT Summary of this function goes here
%   Detailed explanation goes here
    NO_OF_SYMMETRIES = 2;

    
    % TEBD environment initialization:
    EmptyMPS = LSMPS_create(L,NO_OF_SYMMETRIES);
    spinhalf_site = SITE_generate_SU2xU1_spinhalf_fermion_open(SU2sym,U1sym);
    Sites = cellfun(@(x) spinhalf_site,cell(1,L),'UniformOutput',false);  % cell array of spinhalf sites
    TEBD_env = TEBD_init(EmptyMPS,{SU2sym,U1sym},Sites,'FERMIONS',true);

    % Building up the Hamiltonian
    hopp_lr = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'fdag','f',...
                                                 {{[1,2],-1},{[2,1],1}});
    hopp_rl = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'f','fdag',...
                                                 {{[1,2],-1},{[2,1],1}});
    hub_l = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'hub','id',...
                                                 {{[1,1],1}});
    hub_r = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'id','hub',...
                                                 {{[1,1],1}});
    n_l = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'n','id',...
                                                 {{[1,1],1}});


    hopptild_lr = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'ftilddag','ftild',...
                                                 {{[1,2],1},{[2,1],-1}});
    hopptild_rl = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'ftild','ftilddag',...
                                                 {{[1,2],1},{[2,1],-1}});
    hubtild_l = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'hubtild','id',...
                                                 {{[1,1],1}});
    hubtild_r = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'id','hubtild',...
                                                 {{[1,1],1}});
    ntild_l = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'ntild','id',...
                                                 {{[1,1],1}});

    f__ftild_l = COUP_generate_TwoSite_full_OpProd({SU2sym,U1sym},spinhalf_site,spinhalf_site,'f__ftild','id',...
                                                 {{[1,1],1}});
    %hub_r hubtild_r valamiert rossz
    %1st bond
    TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env,1,{{hopp_lr,-1/2},{hopp_rl,-1/2},{hub_l,U},{hub_r,U/2},{hopptild_lr,1/2},{hopptild_rl,1/2},{hubtild_l,-U},{hubtild_r,-U/2},{f__ftild_l,2*G},{n_l,-1i*G},{ntild_l,-1i*G}});
    %last bond
    TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env,L-1,{{hopp_lr,-1/2},{hopp_rl,-1/2},{hub_l,U/2},{hub_r,U},{hopptild_lr,1/2},{hopptild_rl,1/2},{hubtild_l,-U/2},{hubtild_r,-U}});
    for bond_pos = 2:L-2
        TEBD_env = TEBD_set_TwoSite_H_full(TEBD_env,bond_pos,{{hopp_lr,-1/2},{hopp_rl,-1/2},{hub_l,U/2},{hub_r,U/2},{hopptild_lr,1/2},{hopptild_rl,1/2},{hubtild_l,-U/2},{hubtild_r,-U/2}});
    end

    % Determination of two site evolvers (both full and reduced);
    TEBD_env = TEBD_generate_TwoSite_full_evolver_list(TEBD_env,{SU2sym,U1sym},dt);
    %TEBD_env = TEBD_generate_reduced_TwoSite_ops(TEBD_env);
end

