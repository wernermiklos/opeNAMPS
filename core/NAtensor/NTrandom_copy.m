function  out = NTrandom_copy(obj)
            % Returns a random copy of the tensor (same block sizes as in
            % obj, but tensor elements are random (Gaussian) numbers.
            out = NAtensor(obj.leg_names, obj.leg_types, obj.dependencies, obj.no_of_symmetries);
            for key = obj.data.keys
                out.data(key{1}) = randn(size(obj.data(key{1})));
            end
end

