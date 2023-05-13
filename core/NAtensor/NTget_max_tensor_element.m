function maxval = NTget_max_tensor_element(obj)
            % Returns the absolute value of the tensor element with maximal
            % absolute value.
            maxval = 0.0;
            for key = obj.data.keys
                block = obj.data(key{1});
                x = max(abs(block(:)));
                if x > maxval
                    maxval = x;
                end
            end
end

