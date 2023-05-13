function out = LNTget_block_shape(obj, IrrepPositions, IrrepIndices, LegNames)
            % Returns the block's shape
            subshapes = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                subirrepindices = cell(1,length(IrrepIndices));
                for i = 1:length(IrrepIndices)
                    subirrepindices{i} = IrrepIndices{i}(symID);
                end
                subshapes{symID} = NTget_block_shape(obj.sub_tensors{symID},IrrepPositions,subirrepindices,LegNames);
                if isempty(subshapes{symID})
                    out = [];
                    return
                end
            end
            tmp = subshapes{1};
            for symID = 2:obj.no_of_symmetries
                tmp = tmp.*subshapes{symID};
            end
            out = tmp;

end

