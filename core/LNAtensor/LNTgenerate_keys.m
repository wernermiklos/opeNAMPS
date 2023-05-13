function keylist = LNTgenerate_keys(obj, IrrepPositions, IrrepIndices)
            % Generates the key-strings from IrrepPositions and IrrepIndices
            % output is a cell array of length obj.no_of_symmetries. In the
            % cells we put the sub_tensors' keys.
            keylist = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                irrepindices_tmp = cell(1,length(IrrepIndices));
                for i = 1:length(IrrepIndices)
                    irrepindices_tmp{i} = IrrepIndices{i}(symID); 
                end
                keylist{symID} = NTgenerate_key(obj.sub_tensors{symID},IrrepPositions,irrepindices_tmp);
            end


end

