function out = NTlocate_irreps(obj, IrrepPositions)
% Returns the (internal) position of irreps specified in
           % IrrepPositions
           out = zeros(1,length(IrrepPositions));
           for i = 1:length(IrrepPositions)
               legpos = find(strcmp(obj.leg_names,IrrepPositions{i}{1}));
               if isempty(legpos)
                  error(['Error in NAtensor.locate_irreps(): leg "', ...
                      IrrepPositions{i}{1}, '" is not found.'])
               end
               out(i) = obj.dependencies{legpos}(IrrepPositions{i}{2});
           end   
end

