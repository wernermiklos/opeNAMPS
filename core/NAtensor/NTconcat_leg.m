function [out,concat_info] = NTconcat_leg(tensor_list,legname)
%NTCONCAT_LEG Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tensor_list{1}.data)
        error('NTconcat_leg does not work for empty tensors!');
    end
    %First locate the leg in tensor_list{1}
    legpos = NTlocate_legs(tensor_list{1},{legname});
    
    legnum = length(tensor_list{1}.leg_names);
    tenlist_length = length(tensor_list);
    
    no_of_symmetries = tensor_list{1}.no_of_symmetries;

    % We determine the dependencies of the concatted leg, and create the
    % dep_layout (i.e. the set of key elements that the concatted leg depends on.) 
    legdeps = tensor_list{1}.dependencies{legpos};
    dep_layout = zeros(1,length(legdeps)*no_of_symmetries);
    for i = 1:length(legdeps)
        dep_layout(((i-1)*no_of_symmetries+1):(i*no_of_symmetries)) = ((legdeps(i)-1)*no_of_symmetries+1):(legdeps(i)*no_of_symmetries);
    end





    %Create the 'out-map' (key - blockID map), and fill first with the keys
    %of tensor_list{1}
    outmap = containers.Map(keys(tensor_list{1}.data), 1:length(tensor_list{1}.data));
    outkeys = keys(tensor_list{1}.data);



    
    %Create cell array for the outvalues (to concatenate) and fill the
    %first row with the data of tensor_list{1}
    outvalues = cell(length(tensor_list),length(tensor_list{1}.data));
    outvalues(1,:) = values(tensor_list{1}.data);
    
    %Create an array for the dimensions of the concatenated leg, and set
    %its fist element to the dimension in the first tensor (tensor_list{1})
    
    
    %newlegdims = containers.Map('UniformValues',false);
    
    concat_info = containers.Map('UniformValues',false);
    for bID = 1:length(outkeys)
        seckey = outkeys{bID}(dep_layout);
        if ~isKey(concat_info,seckey)
            concat_info(seckey) = [size(outvalues{1,bID},legpos),zeros(1,tenlist_length-1)];
        else
            tmp = concat_info(seckey);
            if tmp(1) ~= size(outvalues{1,bID},legpos)
                error('Dimension mismatch in tensor_list{1}');
            end
        end
    end
    
    
    for i = 2:length(tensor_list)
        if isempty(tensor_list{i}.data)
            error('NTconcat_leg does not work for empty tensors!');
        end
        
        % Get the leg_order and key_layout of tensor_list{i} and transform
        % its keys and blocks accordingly.
        [other_leg_positions, ~, other_key_layout] = NTconsistency_trans(tensor_list{1},tensor_list{i});
        other_keys = cellfun(@(x) x(other_key_layout), keys(tensor_list{i}.data), 'UniformOutput', false);
        other_values = cellfun(@(x) permute(x, other_leg_positions), values(tensor_list{i}.data), 'UniformOutput', false);
        
        %Check which keys are already present in the outmap
        found = isKey(outmap,other_keys);
        for bID_other = 1:length(other_keys)
            seckey = other_keys{bID_other}(dep_layout);
            if isKey(concat_info,seckey)
                concatlegdims_tmp = concat_info(seckey);
            else
                concatlegdims_tmp = zeros(1,tenlist_length);
            end
            if found(bID_other)
                bID = outmap(other_keys{bID_other});
                outvalues{i,bID} = other_values{bID_other};
                concatlegdims_tmp(i) = size(other_values{bID_other},legpos);
                concat_info(seckey) = concatlegdims_tmp;
            else
                bID = size(outvalues,2) + 1;
                outkeys{end+1} = other_keys{bID_other};
                outvalues(:,end+1) = cell(tenlist_length,1);
                outvalues{i,end} = other_values{bID_other};
                %if other_key is not found we must fill the (empty) blocks of
                % outvalues with proper sized zero tensors
                size_tmp = size(outvalues{i,end});
                if legnum > 1
                    size_tmp = [size_tmp,ones(1,legnum-length(size_tmp))];
                end
                concatlegdims_tmp(i) = size_tmp(legpos);
                concat_info(seckey) = concatlegdims_tmp;
                for j = 1:(i-1)
                    size_tmp(legpos) = concatlegdims_tmp(j);
                    outvalues{j,end} = zeros(size_tmp);
                end
            end
        end
        
        %Finally we have to fill the empty blocks in outvalues{i,:} with
        %proper sized zero tensors.
        for bID = 1:length(outkeys)
             seckey = outkeys{bID}(dep_layout);
             concatlegdims_tmp = concat_info(seckey);
             if isempty(outvalues{i,bID})
                size_tmp = size(outvalues{i-1,bID});
                if legnum > 1
                    size_tmp = [size_tmp,ones(1,legnum-length(size_tmp))];
                end
                size_tmp(legpos) = concatlegdims_tmp(i);
                outvalues{i,bID} = zeros(size_tmp);
             end
        end
    end
    
    %We concatenate the tensor blocks and put in the outvalues{1,:} position.
    for bID = 1:size(outvalues,2)
        tmp = outvalues(:,bID);
        outvalues{1,bID} = cat(legpos,tmp{:});
    end
    
    %We create the out tensor and
    out = struct();
    out.type = 'NAtensor';
    out.leg_names = tensor_list{1}.leg_names;
    out.leg_types = tensor_list{1}.leg_types;
    out.no_of_symmetries = tensor_list{1}.no_of_symmetries;
    out.irrep_number = tensor_list{1}.irrep_number;
    out.dependencies = tensor_list{1}.dependencies;
    out.data = containers.Map(outkeys,outvalues(1,:),'UniformValues',false);
    
                    
                
                
            
        



end

