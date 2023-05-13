function out = SUPER_conj_rep(symmetries,Gamma)
%SUPER_CONJ_REP Generates the conj rep for multiple symmetries
%   
out = zeros(1,length(symmetries));
for i = 1:length(symmetries)
    out(i) = symmetries{i}.conj_reps(Gamma(i));
end

end

