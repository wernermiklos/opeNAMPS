function out = LNTrandom_copy(obj)
            % Returns a random copy of the tensor (same subblock sizes as in
            % obj, but tensor elements are random (Gaussian) numbers.
            out = LNAtensor(obj.leg_names, obj.leg_types, obj.dependencies, obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                out = LNTset_subtensor(out,symID,NTrandom_copy(obj.sub_tensors{symID}));
            end


end

