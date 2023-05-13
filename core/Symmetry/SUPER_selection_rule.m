function out = SUPER_selection_rule(symmetries, Gamma1, Gamma2, Gamma)
%SUPER_SELECTION_RULE Summary of this function goes here
%   Detailed explanation goes here
out = 1;
for i = 1:length(symmetries)
    out = out*symmetries{i}.selection_rule(symmetries{i},Gamma1(i),Gamma2(i),Gamma(i));
end
    

end

