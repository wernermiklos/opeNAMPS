function out = LNTget_tensor_element(obj, IrrepPositions, IrrepIndices, LegNames, LegIndices)
            % Returns the tensor element specified by the arguments.
            % IrrepPositions and IrrepIndices specify the block.
            % LegNames and LegIndices specify the tensor element within the
            % block.
            block = LNTget_block(obj,IrrepPositions,IrrepIndices,LegNames);
            shape =  LNTget_block_shape(obj,IrrepPositions,IrrepIndices,LegNames);
            fact = ones(1,length(shape));
            for i = (length(shape)):-1:1
                fact(i) = prod(shape(1:(i-1)));    % some 
            end
            a = reshape(block,[1,prod(shape)]);
            out = a(dot(fact,LegIndices-1)+1);


end

