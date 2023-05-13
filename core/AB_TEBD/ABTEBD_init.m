function ABTEBD_env = ABTEBD_init(startMPS, Symmetries, Sites, varargin)
%ABTEBD_MPS_INIT Summary of this function goes here
%   Detailed explanation goes here
    if ~isempty(startMPS.right_matrices)
        error('startMPS should be left canonical');
    end
    
    parameters = struct();
    parameters.FERMIONS = false;
    parameters.INFINITE = false;
    parameters = parameter_updater(parameters,varargin);

    ABTEBD_env = struct();
    ABTEBD_env.StateMPS = startMPS;
    ABTEBD_env.chain_length = startMPS.chain_length;
    ABTEBD_env.Symmetries = Symmetries;
    ABTEBD_env.no_of_symmetries = length(Symmetries);
    ABTEBD_env.Sites = Sites;
    ABTEBD_env.fermions = parameters.FERMIONS;
    ABTEBD_env.infinite = parameters.INFINITE;
    if ABTEBD_env.infinite
        ABTEBD_env.TwoSiteHlist = cell(1,ABTEBD_env.chain_length);
        ABTEBD_env.Ulist = cell(1,ABTEBD_env.chain_length);
        ABTEBD_env.Uhalflist = cell(1,ABTEBD_env.chain_length);   %to second order Trotter
    else
        ABTEBD_env.TwoSiteHlist = cell(1,ABTEBD_env.chain_length-1);
        ABTEBD_env.Ulist = cell(1,ABTEBD_env.chain_length-1);
        ABTEBD_env.Uhalflist = cell(1,ABTEBD_env.chain_length-1);   %to second order Trotter
    end
    
    
        
        



end

