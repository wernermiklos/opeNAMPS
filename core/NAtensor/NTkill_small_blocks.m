function out = NTkill_small_blocks(obj,threshold)
            % Removes small blocks (max abs element < threshold) from the
            % tensor, and returns the reduced tensor.
            out = NAtensor(obj.leg_names, obj.leg_types, obj.dependencies, obj.no_of_symmetries);
            for key = obj.data.keys
                tmpblock = obj.data(key{1});
                if ~all((abs(tmpblock(:)))< threshold) % we keep then
                    out.data(key{1}) = obj.data(key{1});
                end
            end
end

