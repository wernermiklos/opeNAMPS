function entropy = LSMPS_measure_bond_entropy(obj, symmetries, bondpos)
%LSMPS_MEASURE_BOND_ENTROPY Summary of this function goes here
%   Detailed explanation goes here
    schmidt = obj.schmidt_list{bondpos};
    norm = sqrt(NTdot(schmidt,schmidt,{'t_left','t_right'},{'t_right','t_left'}));
    schmidt = NTmult(schmidt,1/norm);
    [irreps,blocks,~] = NTexport_all_data_blocks(schmidt,{{'t_left',1}},{'t_left','t_right'});
    entropy = 0;
    for bID = 1:length(blocks)
        irrepdim = 1;
        if length(symmetries) ~= length(irreps{bID})
            error('wrong number of symmetries.');
        end
        for symID = 1:length(symmetries)
            irrepdim = irrepdim*symmetries{symID}.irrep_dimensions(irreps{bID}{1}(symID));
        end
        schmidt_values_sq = diag(blocks{bID}).^2;
        entropy = entropy - schmidt_values_sq'*log(schmidt_values_sq) + sum(schmidt_values_sq)*log(irrepdim);
    end

end



