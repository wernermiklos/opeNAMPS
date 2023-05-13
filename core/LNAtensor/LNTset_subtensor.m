function out = LNTset_subtensor(obj,symID,subtensor)
            % Sets obj.sub_tensors{symID} to subtensor. 
            if subtensor.no_of_symmetries ~= 1
                error('Subtensor must have no_of_symmetries=1');
            end
            % We check if subtensor has the proper leg_names, leg_types and dependency_structure: 
            [newsub_leg_positions, ~, newsub_key_layout] = NTconsistency_trans(obj.sub_tensors{symID},subtensor);
            out = LNTcopy(obj);
            out.sub_tensors{symID}.data = containers.Map(); % We empty the old data from sub_tensors{symID}
            %out.sub_tensors{symID}.shape = containers.Map(); % We empty the old data from sub_tensors{symID}
            epsilon = 1.0e-14;     % Small value we will use when decide if every block is 1D 1.0 or not.
            
            
            out.sub_just_ones(symID) = true;
            out.sub_just_nums(symID) = true;
            for keycell = subtensor.data.keys()
                key = keycell{1};
                newkey = key(newsub_key_layout);
                datablock = subtensor.data(key);
                datashape = size(datablock);
                if ~all(datashape == 1)  % there are not just onedimensional blocks.
                    out.sub_just_ones(symID) = false;
                    out.sub_just_nums(symID) = false;
                elseif abs(datablock(1)-1.0) > epsilon % if the blocks are one-dimensional but the tensor element is different from 1.0:
                    out.sub_just_ones(symID) = false;
                end
                out.sub_tensors{symID}.data(newkey) = permute(datablock,newsub_leg_positions);
                %out.sub_tensors{symID}.shape(newkey) = datashape(newsub_leg_positions);
            end            
end

