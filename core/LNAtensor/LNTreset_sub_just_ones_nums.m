function out = LNTreset_sub_just_ones_nums(obj)
            % Resets obj.subs_just_ones and obj.subs_just_nums. Loops over
            % the blocks of the subtensors, and check their shape and
            % values.
            epsilon = 1.0e-15;
            out = LNTcopy(obj);
            out.sub_just_ones = true(1,obj.no_of_symmetries);
            out.sub_just_nums = true(1,obj.no_of_symmetries);
            for symID = 1:out.no_of_symmetries
                for keycell = out.sub_tensors{symID}.data.keys()
                    datablock = out.sub_tensors{symID}.data(keycell{1});
                    datashape = out.sub_tensors{symID}.shape(keycell{1});
                    if ~all(datashape == 1)  % there are not just onedimensional blocks.
                        obj.sub_just_ones(symID) = false;
                        obj.sub_just_nums(symID) = false;
                    elseif abs(datablock(1)-1.0) > epsilon % if the blocks are one-dimensional but the tensor element is different from 1.0:
                        obj.sub_just_ones(symID) = false;
                    end
                end  
            end

end

