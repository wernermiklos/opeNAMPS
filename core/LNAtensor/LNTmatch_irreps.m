function out = LNTmatch_irreps(obj,irrep1,irrep2)
            % Matches two irreps, makes one irrep label of two. (Works only 
            % if the two irreps have equal values in every block.)
            % irrep1 / irrep2 format: {'legname', depID}
            new_sub_tensors = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                new_sub_tensors{symID} = NTmatch_irreps(obj.sub_tensors{symID},irrep1,irrep2);
            end
            out = LNAtensor(obj.leg_names,obj.leg_types,new_sub_tensors{1}.dependencies,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                out = LNTset_subtensor(out,symID,new_sub_tensors{symID});
            end
end

