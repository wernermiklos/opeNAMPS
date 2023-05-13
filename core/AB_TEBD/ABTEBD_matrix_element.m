function out = ABTEBD_matrix_element(ABTEBD_env,BraMPS,site_list,op_names)
%ABTEBD_MEASURE_MATRIX_ELEMENT <other|ops|obj>
%   by: Patrik
%   move to: Patrik_open
%   Detailed explanation goes here

MPS=ABTEBD_env.StateMPS;

for i=1:length(site_list)
    MPS=ABMPS_apply_singlesite_op(MPS,...
                                  site_list{i},...
                                  ABTEBD_env.Sites{site_list{i}}.operators(op_names{i}),...
                                  ABTEBD_env.Sites{site_list{i}}.operator_irreps(op_names{i}),...
                                  ABTEBD_env.Sites,...
                                  ABTEBD_env.Symmetries,...
                                  ABTEBD_env.Sites{site_list{i}}.fermionicities(op_names{i}));
end


out=ABMPS_scalarprod(BraMPS,MPS);
end

