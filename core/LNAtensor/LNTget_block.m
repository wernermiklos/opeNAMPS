function out = LNTget_block(obj, IrrepPositions, IrrepIndices, LegNames)
            % Returns the block 
            subblocks = cell(1,obj.no_of_symmetries);
            subshapes = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                subirrepindices = cell(1,length(IrrepIndices));
                for i = 1:length(IrrepIndices)
                    subirrepindices{i} = IrrepIndices{i}(symID);
                end
                subblocks{symID} = NTget_block(obj.sub_tensors{symID},IrrepPositions,subirrepindices,LegNames);
                subshapes{symID} = NTget_block_shape(obj.sub_tensors{symID},IrrepPositions,subirrepindices,LegNames);
                if isempty(subblocks{symID})
                    out = [];
                    return
                end
            end
            tmp = subblocks{1};
            tmpshape = subshapes{1};
            for symID = 2:obj.no_of_symmetries
                [tmp,tmpshape] = NTblock_kron(tmp,subblocks{symID},tmpshape,subshapes{symID});
            end
            out = tmp;
end

