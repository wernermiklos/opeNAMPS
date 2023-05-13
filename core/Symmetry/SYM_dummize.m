function out = SYM_dummize(obj )
% SYM_dummize : creates a "dummy" symmetry, i.e. deletes Clebsch and
% one_per_dim
% won't work if selection_rule and fusion_rule uses the SYM_****_rule_basic
% function. 
out = obj;
out.type = 'Symmetry DUMMY';
out = rmfield(out,'CGtensor');
out = rmfield(out,'one_per_dim');
end

