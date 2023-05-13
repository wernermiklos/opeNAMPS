function out = SYM_fusion_rule_basic(obj,Gamma1,Gamma2)
%SYM_FUSION_RULE_BASIC Summary of this function goes here
%   Detailed explanation goes here
seclist = NTget_leg_sectors(obj.CGtensor,'alpha');
out = {};
for secID = 1:length(seclist)
    if seclist{secID}{1}{1} == Gamma1 &&  seclist{secID}{1}{2} == Gamma2
        out{end+1} = {seclist{secID}{1}{3},seclist{secID}{2}};
    end
end

end

