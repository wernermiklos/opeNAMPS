function irreps = NTget_active_irrep_values(obj, IrrepPositions)
            % Returns the list of irrep values that encode active blocks.
            irrep_positions = NTlocate_irreps(obj,IrrepPositions);
            if length(irrep_positions) ~= length(unique(irrep_positions)) || ...
                    length(unique(irrep_positions)) ~= obj.irrep_number
                error('Irreps are not clearly specified in IrrepPositions!') 
            end
            key_layout = zeros(obj.no_of_symmetries,obj.irrep_number);
            for i = 1:obj.irrep_number
                key_layout(:,i) = (obj.no_of_symmetries*(irrep_positions(i)-1)+1):(obj.no_of_symmetries*irrep_positions(i));
            end
            irreps = cell(1,obj.data.Count);
            bID = 1;
            for key = obj.data.keys
                irreps{bID} = cell(1,obj.irrep_number);
                for irrepID = 1:obj.irrep_number
                    irreps{bID}{irrepID} = double(key{1}(key_layout(:,irrepID)));
                end
                bID = bID + 1;
            end

end

