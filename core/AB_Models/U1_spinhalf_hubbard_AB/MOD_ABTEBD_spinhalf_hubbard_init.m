function ABTEBD_env = MOD_ABTEBD_spinhalf_hubbard_init(U1Sym, StartMPS, U, dt, varargin)
%MOD_ABTEBD_SPINHALF_HUBBARD_INIT Summary of this function goes here
%   Detailed explanation goes here


    parameters = struct();
    parameters.FERMIONS = true;
    parameters.INFINITE = false;
    parameters = parameter_updater(parameters,varargin);

    spinhalf_fermion_site = MOD_ABSITE_gen_spinhalf_fermion_U1xU1();
    sites = cell(1,StartMPS.chain_length);
    for i = 1:length(sites)
        sites{i} = spinhalf_fermion_site;
    end
    ABTEBD_env = ABTEBD_init(StartMPS,{U1Sym,U1Sym},sites,'FERMIONS',parameters.FERMIONS,'INFINITE',parameters.INFINITE);
    if parameters.INFINITE
        endbond = StartMPS.chain_length;
    else
        endbond = StartMPS.chain_length - 1;
    end
    
    for bondpos = 1:endbond
        if parameters.INFINITE || (bondpos > 1 && bondpos < endbond)
            ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{{'fdag_u', 'f_u', -1}, ...
                                                                 {'f_u', 'fdag_u', 1}, ...
                                                                 {'fdag_d', 'f_d', -1}, ...
                                                                 {'f_d', 'fdag_d', 1}, ...
                                                                 {'id', 'hub', U/2}, ...
                                                                 {'hub', 'id', U/2}});
        elseif bondpos == 1
            ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{{'fdag_u', 'f_u', -1}, ...
                                                                 {'f_u', 'fdag_u', 1}, ...
                                                                 {'fdag_d', 'f_d', -1}, ...
                                                                 {'f_d', 'fdag_d', 1}, ...
                                                                 {'id', 'hub', U/2}, ...
                                                                 {'hub', 'id', U}});
        elseif bondpos == endbond
            ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{{'fdag_u', 'f_u', -1}, ...
                                                                 {'f_u', 'fdag_u', 1}, ...
                                                                 {'fdag_d', 'f_d', -1}, ...
                                                                 {'f_d', 'fdag_d', 1}, ...
                                                                 {'id', 'hub', U}, ...
                                                                 {'hub', 'id', U/2}});
        end
    end
    ABTEBD_env = ABTEBD_gen_Ulist(ABTEBD_env,dt,30);

    
    
end

