function ABTEBD_env = MOD_ABTEBD_U1_nearest_neighbor_spinless_fermion_init(U1sym, StartMPS, V, dt, varargin)
%MOD_ABTEBD_SPINHALF_HUBBARD_INIT Summary of this function goes here
%   
%   H = sum_i [-1*(fdag_{i} f_{i+1} + fdag_{i+1} f_{i}) + V n_{i} n_{i+1}]
%             
%       the sum goes for i = 1:(L-1) if INFINITE is false, and for i = 1:L
%       if INFINITE is true (then L is the length of the unit cell)
%

    parameters = struct();
    parameters.INFINITE = false;
    parameters.Q_EMPTY = -1;
    parameters.Q_OCCUPIED = +1;
    parameters = parameter_updater(parameters,varargin);

    spinless_fermion_site = ABSITE_U1_spinless_fermion(U1sym,parameters.Q_EMPTY, parameters.Q_OCCUPIED);
    sites = cell(1,StartMPS.chain_length);
    for i = 1:length(sites)
        sites{i} = spinless_fermion_site;
    end
    ABTEBD_env = ABTEBD_init(StartMPS,{U1sym},sites,'FERMIONS',true,'INFINITE',parameters.INFINITE);
    if parameters.INFINITE
        endbond = StartMPS.chain_length;
    else
        endbond = StartMPS.chain_length - 1;
    end
    
    for bondpos = 1:endbond
            ABTEBD_env = ABTEBD_set_twositeH(ABTEBD_env,bondpos,{{'fdag', 'f', -1}, ...
                                                                 {'f', 'fdag', 1}, ...     % keep the left-to-right ordering in mind!
                                                                 {'n', 'n', V}});
    end
    ABTEBD_env = ABTEBD_gen_Ulist(ABTEBD_env,dt,30);

    
    
end

