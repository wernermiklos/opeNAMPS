function  out = LNTget_active_sub_irrep_values(obj, IrrepPositions)
            % Returns the list of irrep values (per symmetry) that encode
            % active blocks of subtensors.
            % Output is a cell array of length obj.no_of_symmetries
            out = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                out{symID} = NTget_active_irrep_values(obj.sub_tensors{symID}, IrrepPositions);
            end

end

