function out = LNTcopy(obj)
            % Creates a (deep) copy of the object.
            out = LNAtensor(obj.leg_names, obj.leg_types, obj.dependencies, ...
                           obj.no_of_symmetries);
            out.sub_just_ones = obj.sub_just_ones;
            out.sub_just_nums = obj.sub_just_nums;
            for symID = 1:obj.no_of_symmetries
                out.sub_tensors{symID} = NTcopy(obj.sub_tensors{symID});
            end
end

