function  [newblock,newshape] = NTblock_kron(obj,other,objshape,othershape)
            % Generalization of the kron() function for general (N-dimensional)
            % tensors (arrays).
            %  - obj, other: N-dimensional arrays
            %  - objshape, othershape: the sizes of the arrays (containing also the
            %                          singleton dimensions)
            if length(objshape) ~= length(othershape)
                error('Error: the two tensors must have the same number of legs.')
            end
            newshape = (objshape).*(othershape);
            permute_order = zeros(1,2*length(objshape));
            for i = 1:length(objshape)
                permute_order((2*i-1):(2*i)) = [i, i + length(objshape)];
            end
            newblock = reshape(permute(reshape(obj(:)*(other(:).'),[objshape,othershape]), ...
                                        permute_order),newshape);
end

