function out = SYM_selection_rule_basic(obj,Gamma1,Gamma2,Gamma)
%SYM_SELECTION_RU Summary of this function goes here
%   Detailed explanation goes here
    blockshape = NTget_block_shape(obj.CGtensor,{{'m1',1},{'m2',1},{'M',1}},{Gamma1,Gamma2,Gamma},{'m1','m2','M','alpha'});
    if isempty(blockshape)
        out = 0;
    else
        out = blockshape(4);
    end
end

