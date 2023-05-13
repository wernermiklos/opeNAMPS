function TEBD_env = TEBD_init(startLSMPS, Symmetries, Sites, varargin)
%ABTEBD_MPS_INIT Summary of this function goes here
%   Detailed explanation goes here
    if ~strcmp(startLSMPS.type,'LSMPS')
        error('startMPS should be LeftSymmetric MPS (LSMPS)');
    end
    
    
    parameters = struct();
    % Default parameters
    parameters.FERMIONS = false;
    parameters.INFINITE = false;
    parameters.DIFF_SITES = false;
    parameters.BONDREPS_PER_SYM = {};
    % Update according to varargin
    parameters = parameter_updater(parameters,varargin);


    TEBD_env = struct();
    TEBD_env.StateMPS = startLSMPS;
    TEBD_env.chain_length = startLSMPS.chain_length;
    TEBD_env.Symmetries = Symmetries;
    TEBD_env.no_of_symmetries = length(Symmetries);
    TEBD_env.fermions = parameters.FERMIONS;
    TEBD_env.infinite = parameters.INFINITE;
    TEBD_env.diff_sites = parameters.DIFF_SITES;
    TEBD_env.envelope_list = {};

    if TEBD_env.diff_sites
        error('not implemented yet')
    else
        TEBD_env.Sites = cell(1,TEBD_env.chain_length);
        
        TEBD_env.Sites = cellfun(@(x) Sites{1},TEBD_env.Sites,'UniformOutput',false);    % We make sure that Sites{1} is set to every site.
        envelope = SUPER_generate_TEBD_envelope(Symmetries,Sites{1},Sites{1},...
                                                                 'BONDREPS_PER_SYM', parameters.BONDREPS_PER_SYM);
        TEBD_env.envelope_list = cellfun(@(x) envelope,TEBD_env.Sites, 'UniformOutput',false);  % We make sure that envelope is set to every site.
    end

    if TEBD_env.infinite
        TEBD_env.TwoSiteHlist_full = cell(1,TEBD_env.chain_length);
        TEBD_env.Ulist_full = cell(1,TEBD_env.chain_length);
        TEBD_env.Uhalflist_full = cell(1,TEBD_env.chain_length);   %to second order Trotter
        TEBD_env.TwoSiteHlist_red = cell(1,TEBD_env.chain_length);
        TEBD_env.Ulist_red = cell(1,TEBD_env.chain_length);
        TEBD_env.Uhalflist_red = cell(1,TEBD_env.chain_length);   %to second order Trotter
    else
        TEBD_env.TwoSiteHlist_full = cell(1,TEBD_env.chain_length-1);
        TEBD_env.Ulist_full = cell(1,TEBD_env.chain_length-1);
        TEBD_env.Uhalflist_full = cell(1,TEBD_env.chain_length-1);   %to second order Trotter
        TEBD_env.TwoSiteHlist_red = cell(1,TEBD_env.chain_length-1);
        TEBD_env.Ulist_red = cell(1,TEBD_env.chain_length-1);
        TEBD_env.Uhalflist_red = cell(1,TEBD_env.chain_length-1);   %to second order Trotter
    end
    
    
        
        



end

