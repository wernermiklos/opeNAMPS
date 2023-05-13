function out = ABSITE_anticommutator(obj_in,op1_name,op2_name)
%ABSITE_ANTICOMMUTATOR anticommutator of op1 and op2, {op1,op2}
%   by: Patrik
%   Detailed explanation goes here
out=NTadd(NTdot(obj_in.operators(op1_name),obj_in.operators(op2_name),{'tau'},{'tau~'}),NTdot(obj_in.operators(op2_name),obj_in.operators(op1_name),{'tau'},{'tau~'}));
end

