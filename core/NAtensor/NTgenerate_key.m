function key = NTgenerate_key(obj, IrrepPositions, IrrepIndices)
% Generates the key-string from IrrepPositions and IrrepIndices
            irrep_list = zeros(1,length(IrrepPositions));
            for i = 1:length(IrrepPositions)
                irrep_list(i) = obj.dependencies{strcmp(obj.leg_names, IrrepPositions{i}{1})}(IrrepPositions{i}{2});
            end
            if length(unique(irrep_list)) ~= obj.irrep_number
               error(['Error in NAtensor.generate_key():  IrrepPositions does not', ...
                   ' define all the irrep positions.'])
            end
            ordered_irrepindices = cell(1,length(irrep_list));
            for i = 1:length(irrep_list)
                if length(IrrepIndices{i}) ~=  obj.no_of_symmetries || isa(IrrepIndices,'double')
                   error(['Error in NAtensor.generate_key():  Elements of IrrepIndices', ...
                       'must be an array of obj.no_of_symmetries elements.'])
                end
                ordered_irrepindices{irrep_list(i)} = IrrepIndices{i};
            end
            key = char(zeros(1, obj.no_of_symmetries* length(ordered_irrepindices)));
            for i = 1:length(ordered_irrepindices)
                key((obj.no_of_symmetries*(i-1)+1):obj.no_of_symmetries*i) = char(ordered_irrepindices{i});
            end
end

