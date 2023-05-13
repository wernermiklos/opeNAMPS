function out = LNTset_subblock(obj,symID,IrrepPositions, IrrepIndices, LegNames, datablock)
            % Sets a block of obj.sub_tensors{symID}. It calls the subtensor's set_block function.
             epsilon=1.0e-14;
             out = LNTcopy(obj);
             if ~all(size(datablock) == 1)  % there are not just onedimensional blocks.
                    out.sub_just_ones(symID) = false;
                    out.sub_just_nums(symID) = false;
             elseif abs(datablock(1)-1.0) > epsilon % if the blocks are one-dimensional but the tensor element is different from 1.0:
                    out.sub_just_ones(symID) = false;
             end
             out.sub_tensors{symID} = NTset_block(out.sub_tensors{symID}, IrrepPositions, IrrepIndices, LegNames, datablock);

end

