function Gamma = SUPER_Qnums_to_Gamma(Symmetries,Qnums)
%SUPER_QNUMS Converts Qnums to Gamma for many symmetries.
%   Detailed explanation goes here
    Gamma = zeros(1,length(Symmetries));
    for symID = 1:length(Symmetries)
        Gamma(symID) = Symmetries{symID}.Qnum_to_Gamma(Qnums(symID));
    end
end

