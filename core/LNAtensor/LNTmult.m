function out = LNTmult(obj,val)
            % Returns obj*val. 
            % The function tries to multiply a subtensor that is not
            % "just_ones". If it is not possible, the last subtensor will
            % be multiplied.
            out = LNAtensor(obj.leg_names,obj.leg_types,obj.dependencies,obj.no_of_symmetries);
            out.sub_just_ones = obj.sub_just_ones;
            out.sub_just_nums = obj.sub_just_nums;
            happened = false;
            for symID = 1:obj.no_of_symmetries
                if (~obj.sub_just_ones(symID)) && ~happened && symID < obj.no_of_symmetries
                    out.sub_tensors{symID} = NTmult(obj.sub_tensors{symID},val);
                    happened = true;
                elseif ~happened && symID == obj.no_of_symmetries
                    out.sub_tensors{symID} = NTmult(obj.sub_tensors{symID},val);
                    out.sub_just_ones(symID) = false;
                    happened = true;
                else
                    out.sub_tensors{symID} = NTcopy(obj.sub_tensors{symID});
                end
            end
end

