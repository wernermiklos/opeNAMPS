function varargout = NTslice_leg(obj, legname, concat_info)
%NTSLICE Summary of this function goes here
%   Detailed explanation goes here
    
    no_of_symmetries = obj.no_of_symmetries;
    legpos = NTlocate_legs(obj,{legname});
    legnum = length(obj.leg_names);

    legdeps = obj.dependencies{legpos};
    dep_layout = zeros(1,length(legdeps)*no_of_symmetries);
    for i = 1:length(legdeps)
        dep_layout(((i-1)*no_of_symmetries+1):(i*no_of_symmetries)) = ((legdeps(i)-1)*no_of_symmetries+1):(legdeps(i)*no_of_symmetries);
    end


    %if ~isempty(obj.dependencies{legpos})
    %    error('slicing works only for legs without dependencies')
    %end

    ccinfo_tmp = values(concat_info);
    if isempty(ccinfo_tmp)
        error('concat info should not be empty!');
    end
    no_of_outputs = length(ccinfo_tmp{1});
    keyset = keys(obj.data);
    blocks = values(obj.data);
    newblocks = cell(no_of_outputs,length(blocks));
    for i = 1:length(blocks)
        size_arg = num2cell(size(blocks{i}));  
        seckey = keyset{i}(dep_layout);
        size_arg{legpos} = concat_info(seckey);
        newblocks(:,i) = mat2cell(blocks{i},size_arg{:});
    end

    varargout = cell(1,no_of_outputs);
    for outID = 1:no_of_outputs
        kept = cellfun(@(x) any(x(:)~=0),newblocks(outID,:));
        varargout{outID} = struct();
        varargout{outID}.type = 'NAtensor';
        varargout{outID}.leg_names = obj.leg_names;
        varargout{outID}.leg_types = obj.leg_types;
        varargout{outID}.data = containers.Map(keyset(kept),newblocks(outID,kept),'UniformValues',false);
        varargout{outID}.no_of_symmetries = obj.no_of_symmetries;
        varargout{outID}.irrep_number = obj.irrep_number;
        varargout{outID}.dependencies = obj.dependencies;
    end


end

