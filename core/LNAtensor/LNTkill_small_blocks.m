function out = LNTkill_small_blocks(obj,threshold)
            % Removes small blocks (max abs element < threshold) from the tensor.
            out = LNAtensor(obj.leg_names, obj.leg_types, obj.dependencies, ...
                           obj.no_of_symmetries);
            out.sub_just_ones = obj.sub_just_ones;
            out.sub_just_nums = obj.sub_just_nums;
            
            for symID = 1:obj.no_of_symmetries
                out = LNTset_subtensor(out,symID,NTkill_small_blocks(obj.sub_tensors{symID},threshold));
            end

end

