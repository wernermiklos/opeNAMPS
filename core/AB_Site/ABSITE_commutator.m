function out = ABSITE_commutator(obj,op1_name,op2_name)
%ABSITE_COMMUTATOR commutator of op1 and op2, [op1,op2]
%   by: Patrik
%   Detailed explanation goes here
out=NTsubtr(NTdot(obj.operators(op1_name),obj.operators(op2_name),{'tau'},{'tau~'}),NTdot(obj.operators(op2_name),obj.operators(op1_name),{'tau'},{'tau~'}));
end
