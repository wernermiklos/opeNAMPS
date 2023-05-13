function out = SUPER_fusion_rule(symmetries,Gamma1,Gamma2 )
%SUPER_FUSION_RULE Summary of this function goes here
%   Detailed explanation goes here
    tmp_old = symmetries{1}.fusion_rule(symmetries{1},Gamma1(1),Gamma2(1));
    for i = 2:length(symmetries)
        fr_tmp = symmetries{i}.fusion_rule(symmetries{i},Gamma1(i),Gamma2(i));
        tmp = cell(1,length(tmp_old)*length(fr_tmp));
        j = 1;
        for secID1 = 1:length(tmp_old)
            for secID2 = 1:length(fr_tmp)
                tmp{j} = {[tmp_old{secID1}{1},fr_tmp{secID2}{1}], tmp_old{secID1}{2}* fr_tmp{secID2}{2}};
                j = j + 1;
            end
        end
        tmp_old = tmp;
    end
    out = tmp_old;
end
