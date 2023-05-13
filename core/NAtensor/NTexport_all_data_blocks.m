function [irreps, blocks, shapes] = NTexport_all_data_blocks(obj, IrrepPositions, LegNames)
            % Exports list of all data blocks and the corresponding irrep values in cell arrays.
            % The irrep positions are specified in IrrepPositions
            % {{'legname1',depID1},{'legname2',depID2},...}, the desired
            % leg order is specified in LegNames.
            irrep_positions = NTlocate_irreps(obj,IrrepPositions);
            if length(irrep_positions) ~= length(unique(irrep_positions)) || ...
                    length(unique(irrep_positions)) ~= obj.irrep_number
                error('Irreps are not clearly specified in IrrepPositions!') 
            end
            leg_positions = NTlocate_legs(obj,LegNames);
            if length(leg_positions) ~= length(unique(leg_positions)) || ...
                    length(unique(leg_positions)) ~= length(obj.leg_names)
                 error('Legs are not clearly specified in IrrepPositions!')
            end
            key_layout = zeros(obj.no_of_symmetries,obj.irrep_number);
            for i = 1:obj.irrep_number
                key_layout(:,i) = (obj.no_of_symmetries*(irrep_positions(i)-1)+1):(obj.no_of_symmetries*irrep_positions(i));
            end
            irreps = cell(1,obj.data.Count);
            blocks = cell(1,obj.data.Count);
            shapes = cell(1,obj.data.Count);
            legnum = length(obj.leg_names);
            bID = 1;
            for key = obj.data.keys
                blocks{bID} = permute(obj.data(key{1}),leg_positions);
                tmpshape = [size(blocks{bID}),ones(1,legnum-ndims(blocks{bID}))];
                if legnum == 1
                    shapes{bID} = tmpshape(1);
                else
                    shapes{bID} = tmpshape;   %tmpshape(leg_positions);   WAS WRONG!!!
                end
                irreps{bID} = cell(1,obj.irrep_number);
                for irrepID = 1:obj.irrep_number
                    irreps{bID}{irrepID} = double(key{1}(key_layout(:,irrepID)));
                end
                bID = bID + 1;
            end


end

