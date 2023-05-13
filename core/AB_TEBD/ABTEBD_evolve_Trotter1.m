function [ABTEBD_env, info] = ABTEBD_evolve_Trotter1(ABTEBD_env, Mmax, Nstep,varargin)
%ABTEBD_EVOLVE_TROTTER1 Summary of this function goes here
%   Detailed explanation goes here

parameters = struct();
parameters.NORMALIZE = true;
parameters = parameter_updater(parameters,varargin);

info = struct();
info.bond_dims = zeros(1,length(ABTEBD_env.Ulist));
info.trunc_errs = zeros(1,length(ABTEBD_env.Ulist));
for i = 1:Nstep
    for pos = 1:2:length(ABTEBD_env.Ulist)
        [ABTEBD_env, step_info] = ABTEBD_evolve2site(ABTEBD_env,pos,ABTEBD_env.Ulist{pos}, Mmax, 'NORMALIZE', parameters.NORMALIZE);
        info.bond_dims(pos) = step_info.bond_dim;
        info.trunc_errs(pos) = step_info.trunc_err;
    end
    for pos = 2:2:length(ABTEBD_env.Ulist)
        [ABTEBD_env, step_info] = ABTEBD_evolve2site(ABTEBD_env,pos,ABTEBD_env.Ulist{pos}, Mmax, 'NORMALIZE', parameters.NORMALIZE);
        info.bond_dims(pos) = step_info.bond_dim;
        info.trunc_errs(pos) = step_info.trunc_err;
    end
end
    
    
end

