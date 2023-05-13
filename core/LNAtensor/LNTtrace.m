function out = LNTtrace(obj,LegNames1,LegNames2)
            % Calculates the trace of the tensor, i.e. contracts legs
            % listed in LegNames1 with the ones listed in LegNames2
            new_sub_tensors = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                new_sub_tensors{symID} = NTtrace(obj.sub_tensors{symID},LegNames1,LegNames2);
            end
            out = LNAtensor(new_sub_tensors{1}.leg_names, ...
                            new_sub_tensors{1}.leg_types, ...
                            new_sub_tensors{1}.dependencies, ...
                            obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                out = NTset_subtensor(out, symID,new_sub_tensors{symID});
            end


end

