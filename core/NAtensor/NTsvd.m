function [left_tensor,schmidt_tensor,right_tensor,norm_sq, M] = NTsvd(obj,left_leg_names,right_leg_names,bond_leg_names,bond_leg_types, Mmax,eps)
            % Factorizes the tensor by singular value decomposition. Results in three tensors: "left_tensor",
            % "schmidt_tensor", and "right_tensor". The "left_tensor" contains all the legs specified in LeftLegs, while the
            % "right_tensor" contains legs in RightLegs. (All legs of "self" should be either in LeftLegs or RightLegs). In
            % addition, the "left_tensor" and "right_tensor" have one-one additional "bond-legs" (name specified in
            % BondLegNames) that connect them to the "schmidt_tensor". The number of irreps in "left_tensor", "right_tensor"
            % and "schmidt_tensor" depends on the dependency structure of "self".
            % ----
            % left_leg_names:     List of leg names that will be in the "left_tensor".
            % right_leg_names:    List of leg names that will be in the "right_tensor".
            %                   REMARK: LeftLegNames and RightLegNames should contain together all the legs of "self".
            %                   Otherwise error is raised.
            % bond_leg_names:     Names of bond legs. Format: ['left_bond_name', 'right_bond_name']. "left_tensor" will have a
            %                   leg with name 'left_bond_name', while
            %                   "right_tensor" will have one with
            %                   'right_bond_name'.in
            %                   The (diagonal) Schmidt tensor will have two legs, their names are 'left_bond_name' and
            %                   'right_bond_name'. 
            if nargin == 5
                Mmax = inf;
                eps = 10^(-16);
            elseif nargin == 6
                 eps = 10^(-16);
            elseif nargin > 7
                error('nargin is 4, 5 or 6. (5th optional arg is the max bond dimension Mmax. 6th optional arg is eps (minimal value of kept Schmidt values).');
            end
            if length(unique([left_leg_names,right_leg_names])) ~= length(obj.leg_names) || ...
                    length(unique([left_leg_names,right_leg_names])) ~= length([left_leg_names,right_leg_names])
                error('Error in left_leg_names or right_leg_names.') 
            end
            if length(unique([left_leg_names, bond_leg_names{1}])) ~= length(left_leg_names) + 1 || ...
                length(unique([right_leg_names, bond_leg_names{2}])) ~= length(right_leg_names) + 1 || ...
                strcmp(bond_leg_names{1},bond_leg_names{2})
                error('Error: the new tensors have ambiguously defined leg names')
            end
            left_leg_pos = NTlocate_legs(obj,left_leg_names);
            right_leg_pos = NTlocate_legs(obj,right_leg_names);
            if strcmp(bond_leg_types{1},'i')
                left_leg_types = [obj.leg_types(left_leg_pos), 'o'];  % The last leg of the left tensor is the new bond
            elseif strcmp(bond_leg_types{1},'o')
                left_leg_types = [obj.leg_types(left_leg_pos), 'i'];  % The last leg of the left tensor is the new bond
            else
                error('Leg types are i or o.')
            end
            if strcmp(bond_leg_types{2},'i')
                right_leg_types = ['o',obj.leg_types(right_leg_pos)]; % The first leg of the right tensor is the new bond
            elseif strcmp(bond_leg_types{2},'o')    
                right_leg_types = ['i',obj.leg_types(right_leg_pos)]; % The first leg of the right tensor is the new bond
            else
                error('Leg types are i or o.')
            end
            new_left_leg_names = [left_leg_names, bond_leg_names{1}]; % The last leg of the left tensor is the new bond 
            new_right_leg_names = [bond_leg_names{2}, right_leg_names];  % The first leg of the right tensor is the new bond
            
            left_irrepIDs = sort(unique(cell2mat(obj.dependencies(left_leg_pos))));
            right_irrepIDs = sort(unique(cell2mat(obj.dependencies(right_leg_pos))));
            center_irrepIDs = sort(intersect(left_irrepIDs,right_irrepIDs));
            
            left_key_layout = zeros(1,obj.no_of_symmetries*length(left_irrepIDs));
            for i = 1:length(left_irrepIDs)
                left_key_layout((obj.no_of_symmetries*(i-1)+1):(obj.no_of_symmetries*i)) = ...
                    (obj.no_of_symmetries*(left_irrepIDs(i)-1)+1):(obj.no_of_symmetries*left_irrepIDs(i));
            end
            
            right_key_layout = zeros(1,obj.no_of_symmetries*length(right_irrepIDs));
            for i = 1:length(right_irrepIDs)
                right_key_layout((obj.no_of_symmetries*(i-1)+1):(obj.no_of_symmetries*i)) = ...
                    (obj.no_of_symmetries*(right_irrepIDs(i)-1)+1):(obj.no_of_symmetries*right_irrepIDs(i));
            end
            
            center_key_layout = zeros(1,obj.no_of_symmetries*length(center_irrepIDs));
            for i = 1:length(center_irrepIDs)
                center_key_layout((obj.no_of_symmetries*(i-1)+1):(obj.no_of_symmetries*i)) = ...
                    (obj.no_of_symmetries*(center_irrepIDs(i)-1)+1):(obj.no_of_symmetries*center_irrepIDs(i));
            end
            
            left_dependencies = cell(1,length(left_leg_names)+1);
            for i = 1:length(left_leg_pos)
                legID = left_leg_pos(i);
                left_dependencies{i} = zeros(1,length(obj.dependencies{legID}));
                for depID = 1:length(obj.dependencies{legID})
                    left_dependencies{i}(depID) = find(left_irrepIDs == obj.dependencies{legID}(depID));
                end
            end
            left_dependencies{end} = zeros(1,length(center_irrepIDs));
            for depID = 1:length(center_irrepIDs)
                left_dependencies{end}(depID) = find(left_irrepIDs == center_irrepIDs(depID));
            end
            right_dependencies = cell(1,length(right_leg_names)+1);
            for i = 2:(length(right_leg_pos)+1)
                legID = right_leg_pos(i-1);
                right_dependencies{i} = zeros(1,length(obj.dependencies{legID}));
                for depID = 1:length(obj.dependencies{legID})
                    right_dependencies{i}(depID) = find(right_irrepIDs == obj.dependencies{legID}(depID));
                end
            end
            right_dependencies{1} = zeros(1,length(center_irrepIDs));
            for depID = 1:length(center_irrepIDs)
                right_dependencies{1}(depID) = find(right_irrepIDs == center_irrepIDs(depID));
            end
            center_dependencies = {1:length(center_irrepIDs),1:length(center_irrepIDs)};
            
            left_tensor = NAtensor(new_left_leg_names, left_leg_types, left_dependencies, obj.no_of_symmetries);
            right_tensor = NAtensor(new_right_leg_names, right_leg_types, right_dependencies, obj.no_of_symmetries);
            schmidt_tensor = NAtensor(bond_leg_names, bond_leg_types, center_dependencies, obj.no_of_symmetries);
            
            
            svd_block_centerkeys = {};
            keys_in_svd_blocks = {};
            left_sector_keys = {};
            right_sector_keys = {};
            left_sector_dims = {};
            right_sector_dims = {};
            left_sector_positions = {};
            right_sector_positions = {};
            left_sector_shapes = {};
            right_sector_shapes = {};
            left_totdim = {};
            right_totdim = {};
            objlegnum = length(obj.leg_names);
            
            
            
            for keycell = obj.data.keys
                key = keycell{1};
                leftkey = key(left_key_layout);
                rightkey = key(right_key_layout);
                centerkey = key(center_key_layout);
                datablock = obj.data(key);
                datashape = [size(datablock),ones(1,objlegnum-ndims(datablock))];
                left_sector_shape = datashape(left_leg_pos);
                right_sector_shape = datashape(right_leg_pos);
                left_sector_dim = prod(left_sector_shape);
                right_sector_dim = prod(right_sector_shape);
                if ~any(strcmp(svd_block_centerkeys,centerkey))   % If centerkey appears first
                    svd_block_centerkeys{end + 1} = centerkey;
                    svd_bID = length(svd_block_centerkeys);
                    keys_in_svd_blocks{end + 1} = {};
                    left_sector_keys{end + 1} = {};
                    right_sector_keys{end + 1} = {};
                    left_sector_dims{end + 1} = [];
                    right_sector_dims{end + 1} = [];
                    left_sector_positions{end + 1} = {};
                    right_sector_positions{end + 1} = {};
                    left_sector_shapes{end + 1} = {};
                    right_sector_shapes{end + 1} = {};
                    left_totdim{end + 1} = 0;
                    right_totdim{end + 1} = 0;
                else
                    svd_bID = find(strcmp(svd_block_centerkeys,centerkey));
                end
                
                keys_in_svd_blocks{svd_bID}{end+1} = key;
                
                if ~any(strcmp(left_sector_keys{svd_bID},leftkey))   % if leftkey appears first in the svd_block
                    left_sector_keys{svd_bID}{end + 1} = leftkey;
                    left_secID = length(left_sector_keys{svd_bID});
                    left_sector_dims{svd_bID}(end + 1) = left_sector_dim;
                    left_sector_positions{svd_bID}{end + 1} = (left_totdim{svd_bID}+1):(left_totdim{svd_bID} + left_sector_dim);
                    left_sector_shapes{svd_bID}{end + 1} = left_sector_shape;
                    left_totdim{svd_bID} = left_totdim{svd_bID} + left_sector_dim;
                else
                    left_secID = find(strcmp(left_sector_keys{svd_bID},leftkey));
                    if ~isequal(left_sector_shapes{svd_bID}{left_secID},left_sector_shape)
                        error('NAtensor is corrupted!')
                    end
                end
                
                if ~any(strcmp(right_sector_keys{svd_bID},rightkey))   % if leftkey appears first in the svd_block
                    right_sector_keys{svd_bID}{end + 1} = rightkey;
                    right_secID = length(right_sector_keys{svd_bID});
                    right_sector_dims{svd_bID}(end + 1) = right_sector_dim;
                    right_sector_positions{svd_bID}{end + 1} = (right_totdim{svd_bID}+1):(right_totdim{svd_bID} + right_sector_dim);
                    right_sector_shapes{svd_bID}{end + 1} = right_sector_shape;
                    right_totdim{svd_bID} = right_totdim{svd_bID} + right_sector_dim;
                else
                    right_secID = find(strcmp(right_sector_keys{svd_bID},rightkey));
                    if ~isequal(right_sector_shapes{svd_bID}{right_secID},right_sector_shape)
                        error('NAtensor is corrupted!')
                    end
                end
            end
            
            if Mmax < inf
                schmidt_store = [];
                svd_bID_store = [];
            end
            left_truncated = containers.Map();
            right_truncated = containers.Map();
            
            
            for svd_bID = 1:length(svd_block_centerkeys)
                svd_block = zeros(left_totdim{svd_bID},right_totdim{svd_bID});
                centerkey = svd_block_centerkeys{svd_bID};
                for key_cell = keys_in_svd_blocks{svd_bID}
                    key = key_cell{1};
                    leftkey = key(left_key_layout);
                    rightkey = key(right_key_layout);
                    left_truncated(leftkey) = false;          % Necessary, if Mmax < inf, i.e. if we truncate
                    right_truncated(rightkey) = false;        % Necessary, if Mmax < inf, i.e. if we truncate
                    left_secID = find(strcmp(left_sector_keys{svd_bID},leftkey));
                    right_secID = find(strcmp(right_sector_keys{svd_bID},rightkey));
                    datablock = obj.data(key);
                    matrixblock = reshape(permute(datablock,[left_leg_pos,right_leg_pos]), ...
                        [left_sector_dims{svd_bID}(left_secID),right_sector_dims{svd_bID}(right_secID)]);
                    svd_block(left_sector_positions{svd_bID}{left_secID}, right_sector_positions{svd_bID}{right_secID}) = matrixblock;   
                end
                [U, S, V] = svd(svd_block,'econ');
                V = V';
                schmidt_tensor.data(centerkey) = S;
                schmidt_dim = size(S);
                schmidt_dim = schmidt_dim(1);
                if Mmax < inf
                    schmidt_store = [schmidt_store,diag(S)'];
                    svd_bID_store = [svd_bID_store,svd_bID*ones(1,schmidt_dim)];
                end
                for left_secID = 1:length(left_sector_keys{svd_bID})
                    leftkey = left_sector_keys{svd_bID}{left_secID};
                    left_tensor.data(leftkey) = reshape(U(left_sector_positions{svd_bID}{left_secID},:), ...
                                                               [left_sector_shapes{svd_bID}{left_secID}, schmidt_dim]);
                    
                end
                for right_secID = 1:length(right_sector_keys{svd_bID})
                    rightkey = right_sector_keys{svd_bID}{right_secID};
                    right_tensor.data(rightkey) = reshape(V(:,right_sector_positions{svd_bID}{right_secID}), ...
                                                               [schmidt_dim,right_sector_shapes{svd_bID}{right_secID}]);
                    
                end
                
            end
            leftlegnum = length(left_tensor.leg_names);
            rightlegnum = length(right_tensor.leg_names);
            if Mmax < inf
                Mmax = min(Mmax,sum(schmidt_store >= sqrt(eps)));
                %disp([length(schmidt_store),Mmax]);
                if Mmax < length(schmidt_store)
                    [~,sortIDx] = sort(schmidt_store,'descend');
                    svd_bID_store = svd_bID_store(sortIDx);
                    svd_bID_store = svd_bID_store(1:Mmax);
                    for svd_bID = 1:length(svd_block_centerkeys)
                        centerkey = svd_block_centerkeys{svd_bID};
                        new_schmidt_dim = sum(svd_bID_store == svd_bID);
                        schmidt_block = schmidt_tensor.data(centerkey);
                        new_schmidt_block = schmidt_block(1:new_schmidt_dim,1:new_schmidt_dim);
                        if new_schmidt_dim > 0
                            schmidt_tensor.data(centerkey) = new_schmidt_block;
                        else
                            remove(schmidt_tensor.data,centerkey);
                        end
                        for key_cell = keys_in_svd_blocks{svd_bID}
                            key = key_cell{1};
                            leftkey = key(left_key_layout);
                            rightkey = key(right_key_layout);
                            
                            if left_truncated(leftkey) == false
                                left_truncated(leftkey) = true;
                                left_block = left_tensor.data(leftkey);
                                left_shape = [size(left_block),ones(1,leftlegnum-ndims(left_block))];
                                tmp = reshape(left_block,[prod(left_shape(1:end-1)),left_shape(end)]);
                                tmp = tmp(:,1:new_schmidt_dim);
                                new_left_block = reshape(tmp,[left_shape(1:end-1),new_schmidt_dim]);
                                if new_schmidt_dim > 0
                                    left_tensor.data(leftkey) = new_left_block;
                                else
                                    remove(left_tensor.data,leftkey);
                                end
                            end
                            
                            if right_truncated(rightkey) == false
                                right_truncated(rightkey) = true;
                                right_block = right_tensor.data(rightkey);
                                right_shape = [size(right_block),ones(1,rightlegnum-ndims(right_block))];
                                tmp = reshape(right_block,[right_shape(1),prod(right_shape(2:end))]);
                                tmp = tmp(1:new_schmidt_dim,:);
                                new_right_block = reshape(tmp,[new_schmidt_dim,right_shape(2:end)]);
                                if new_schmidt_dim > 0
                                    right_tensor.data(rightkey) = new_right_block;
                                else
                                    remove(right_tensor.data,rightkey);
                                end
                            end
                        end
                    end
                end
            end
            norm_sq = NTdot(schmidt_tensor,NTconj(schmidt_tensor),bond_leg_names,bond_leg_names);
            if Mmax < inf
                M = min(min(Mmax,sum(schmidt_store >= sqrt(eps))));
            else
                M = schmidt_dim;
            end
end

