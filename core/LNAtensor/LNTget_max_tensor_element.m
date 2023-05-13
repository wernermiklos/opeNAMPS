function maxval = LNTget_max_tensor_element(obj)
            % Returns the absolute value of the tensor element with maximal
            % absolute value.
            maxval = 1.0;
            for symID = 1:obj.no_of_symmetries
                maxval = maxval*NTget_max_tensor_element(obj.sub_tensors{symID}); 
            end
end


