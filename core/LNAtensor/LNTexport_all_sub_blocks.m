function out = LNTexport_all_sub_blocks(obj, IrrepPositions, LegNames)
            % Calls export_all_data_blocks for each subtensors, the result
            % is put into a cell array. 
            out = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                [irreps,blocks] = NTexport_all_data_blocks(obj.sub_tensors{symID}, IrrepPositions,LegNames);
                out{symID} = {irreps,blocks};
            end
end

