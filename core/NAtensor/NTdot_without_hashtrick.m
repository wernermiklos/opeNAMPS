function out = NTdot(obj, other, legnames1, legnames2, renaming)
            % Calculates dot product (leg contraction) of "self" and "other". The contracted legs are specified in
            % "LegNames1" and "LegNames2". In Renaming1 and Renaming2 one can specify new names of (remaining) legs after
            % the contraction. If a leg is not listed in Renaming1 or Renaming2, the old name will be used.
            % If the list LegNames1 and LegNames2 is empty, NAtensor.Simplemult() will be called (Renamings are performed!)
            % ----
            % other:            the other NAtensor.
            % legnames1         list of leg names of "obj" that are contracted.
            % legnames2         list of leg names of "other" that are contracted with the legs of "obj". Length of
            %                   legnames2 must equal length of legnames1. Contractions are performed pairwise (i.e.
            %                   legnames1{0} is contracted with legnames2{0}, ...). The contracted legs must be of the same
            %                   kind, (same dependency structure, same dimensions in the two tensors).
            % renaming          List of leg renamings in "obj" and "other. 
            %                   renaming{1} contains the renamings of "obj"
            %                   renaming{2} containes the renamings of
            %                   "other"
            %                   Format: {{{'oldname1_1', 'newname1_1'}, ...
            %                             {'oldname1_2', 'newname1_2'}, ...},
            %                            {{'oldname2_1', 'newname2_1'}, ...
            %                             {'oldname2_2', 'newname2_2'}}}
            if nargin == 4
                renaming = {{},{}};
            elseif nargin == 5
            else
                error('NTdot() has 3 or 4 parameters. The 4th (optional) parameter is "renaming"')
            end
            if strcmp(other.type,'LNAtensor') || strcmp(obj.type,'LNAtensor')   % to call the LNAtensor.dot() function
                out = LNTdot(obj,other,legnames1,legnames2,renaming);
                return
            end
                
         
            if ~strcmp(other.type,'NAtensor')
                error('Error in NTdot(): other must be an NAtensor too.')
            end
            
            if length(legnames1) ~= length(legnames2)
                error('The number of multiplied legs are not the same in the two tensor!')
            end
            if obj.no_of_symmetries ~= other.no_of_symmetries
                error('The number of symmetries of the two tensors are not equal')
            end
            legpos1 = NTlocate_legs(obj, legnames1);
            legpos2 = NTlocate_legs(other, legnames2);
            freelegs1 = 1:length(obj.leg_names);
            freelegs1(legpos1) = [];
            freelegs2 = 1:length(other.leg_names);
            freelegs2(legpos2) = [];
            
            if length(unique(legpos1)) ~= length(legpos1) || ...
                   length(unique(legpos2)) ~= length(legpos2)
                error('Repetition in legnames1 or legnames2')
            end
            
            for i = 1:length(legnames1)
                if obj.leg_types{legpos1(i)} == other.leg_types{legpos2(i)}
                    error('Legs with different types ("i","o") can only be contracted.')
                end
            end
            keptlegs1 = 1:length(obj.leg_names);
            keptlegs1(legpos1) = [];
            keptlegs2 = 1:length(other.leg_names);
            keptlegs2(legpos2) = [];
            
            % We match the irreps according to the dependencies
            tmpirrepindices1 = 1:obj.irrep_number;
            tmpirrepindices2 = (1:other.irrep_number) + obj.irrep_number;
            matchkey_layout1 = [];
            matchkey_layout2 = [];
            for i = 1:length(legnames1)
                for depID = 1:length(obj.dependencies{legpos1(i)})
                    irrepind1 = obj.dependencies{legpos1(i)}(depID);
                    matchkey_layout1 = [matchkey_layout1, (obj.no_of_symmetries*(irrepind1-1) +1):(obj.no_of_symmetries*irrepind1)];
                    irrepind2 = other.dependencies{legpos2(i)}(depID);
                    matchkey_layout2 = [matchkey_layout2, (other.no_of_symmetries*(irrepind2-1) +1):(other.no_of_symmetries*irrepind2)];
                    if tmpirrepindices1(irrepind1) < tmpirrepindices2(irrepind2)
                        % We change all the index values newirrepindices2(irrepind2) to newirrepindices1(irrpeind1)
                        tmpirrepindices1(tmpirrepindices1 == tmpirrepindices2(irrepind2)) = tmpirrepindices1(irrepind1);
                        tmpirrepindices2(tmpirrepindices2 == tmpirrepindices2(irrepind2)) = tmpirrepindices1(irrepind1);
                    elseif tmpirrepindices1(irrepind1) > tmpirrepindices2(irrepind2)
                        tmpirrepindices2(tmpirrepindices2 == tmpirrepindices1(irrepind1)) = tmpirrepindices2(irrepind2);
                        tmpirrepindices1(tmpirrepindices1 == tmpirrepindices1(irrepind1)) = tmpirrepindices2(irrepind2);
                    end
                end     
            end
            tmpirrepinds_all = unique([tmpirrepindices1,tmpirrepindices2]);
            
            % We sort the irrepindices to make the kept ones first and the
            % summed up after
            
            kept_tmpirrepinds = [];
            for i = 1:length(keptlegs1)
                kept_tmpirrepinds = union(kept_tmpirrepinds,tmpirrepindices1(obj.dependencies{keptlegs1(i)}),'sorted');
            end
            for i = 1:length(keptlegs2)
                kept_tmpirrepinds = union(kept_tmpirrepinds,tmpirrepindices2(other.dependencies{keptlegs2(i)}),'sorted');
            end
            sumup_tmpirrepinds = setdiff(tmpirrepinds_all,kept_tmpirrepinds,'sorted');
            neworder = [kept_tmpirrepinds,sumup_tmpirrepinds];
            new_irrep_number = length(kept_tmpirrepinds);   % the number of irreps in the product tensor
            
            newirrepindices1 = tmpirrepindices1;
            newirrepindices2 = tmpirrepindices2;
            for i = 1:length(neworder)
                newirrepindices1(tmpirrepindices1 == neworder(i)) = i;
                newirrepindices2(tmpirrepindices2 == neworder(i)) = i;
            end
            
            % We may want to do this
%             if length(unique(newirrepindices1)) ~= length(newirrepindices1)    % there is a irrep match within the tensor 'obj' ==> we can 'pre-filter' the blocks
%                 prefilter1 = true;
%             else
%                 prefilter1 = false;
%             end
%             
%             if length(unique(newirrepindices2)) ~= length(newirrepindices2)    % there is a irrep match within the tensor 'other' ==> we can 'pre-filter' the blocks
%                 prefilter2 = true;
%             else
%                 prefilter2 = false;
%             end
            
            % Now we create the list of new leg names. If the name of a free (remaining) leg is listed in renaming{1} (legs
            % of obj), or renaming{2} (legs of other), then the new name is stored. Otherwise we store the name of the leg as
            % it is used in tensors obj or other. If there is ambiguity in the new leg names (the same name is defined for
            % two or more legs) error will be raised later.
            % At the same time we also create the list of new leg types.
            old_leg_names1 = obj.leg_names(keptlegs1);
            new_leg_names1 = old_leg_names1;
            for rnm = renaming{1}
                new_leg_names1{strcmp(rnm{1}{1},old_leg_names1)} = rnm{1}{2};
            end
            
            old_leg_names2 = other.leg_names(keptlegs2);
            new_leg_names2 = old_leg_names2;
            for rnm = renaming{2}
                new_leg_names2{strcmp(rnm{1}{1},old_leg_names2)} = rnm{1}{2};
            end
            
            new_leg_names = [new_leg_names1, new_leg_names2];
            new_leg_types = [obj.leg_types(keptlegs1), other.leg_types(keptlegs2)];
            
            if isempty(new_leg_names)
                just_a_number = true;
            else
                just_a_number = false;
            end
            
            % The dependencies of the new tensor are now determined. 
            if ~just_a_number
                new_dependencies1 = cell(1,length(keptlegs1));
                for i = 1:length(keptlegs1)
                    new_dependencies1{i} = newirrepindices1(obj.dependencies{keptlegs1(i)});
                end
                new_dependencies2 = cell(1,length(keptlegs2));
                for i = 1:length(keptlegs2)
                    new_dependencies2{i} = newirrepindices2(other.dependencies{keptlegs2(i)});
                end
                new_dependencies = [new_dependencies1, new_dependencies2];
            end
            
            
            % We determine how the new block key is generated from key1 and
            % key2
            newkey_layout = [];
            tmp_newirrepindices = [newirrepindices1, newirrepindices2];
            for newirrepID = 1:new_irrep_number
                poslist = find(tmp_newirrepindices == newirrepID);
                pos = poslist(1);
                newkey_layout = [newkey_layout,(obj.no_of_symmetries*(pos-1)+1):(obj.no_of_symmetries*pos)];
            end
                    
            % we create the out variable
            if ~just_a_number
                out = NAtensor(new_leg_names, new_leg_types, new_dependencies, obj.no_of_symmetries);
            else
                out = 0.0;
            end
            
            %permute_order1 = [keptlegs1, legpos1];
            %permute_order2 = [legpos2, keptlegs2];
            
            matchcomb_ids = containers.Map('KeyType','char','ValueType','double');
            
            %transformed_blocks = {};
            %kept_shapes = {};
            mID = 1;
            objkeys = obj.data.keys;
            objvalues = obj.data.values;
            otherkeys = other.data.keys;
            othervalues = other.data.values;
            
            if length(objkeys) < length(otherkeys)
                match_combinations = cell(length(objkeys),2);
                for bID1 = 1:length(objkeys)
                    key_tmp = objkeys{bID1}(matchkey_layout1);
                    if isKey(matchcomb_ids,key_tmp)
                            match_combinations{matchcomb_ids(key_tmp),1}(end+1) = bID1;
                    else
                        matchcomb_ids(key_tmp) = mID;
                        mID = mID + 1;
                        match_combinations{mID-1,1} = bID1;
                    end      
                end     
                for bID2 = 1:length(otherkeys)
                    key_tmp = otherkeys{bID2}(matchkey_layout2);
                    if isKey(matchcomb_ids,key_tmp)
                        match_combinations{matchcomb_ids(key_tmp),2}(end+1) = bID2;
                    end        
                end
            else    
                match_combinations = cell(length(otherkeys),2);
                for bID2 = 1:length(otherkeys)
                    key_tmp = otherkeys{bID2}(matchkey_layout2);
                    if isKey(matchcomb_ids,key_tmp)
                            match_combinations{matchcomb_ids(key_tmp),2}(end+1) = bID2;
                    else
                        matchcomb_ids(key_tmp) = mID;
                        mID = mID + 1;
                        match_combinations{mID-1,2} = bID2;
                    end      
                end
                for bID1 = 1:length(objkeys)
                    key_tmp = objkeys{bID1}(matchkey_layout1);
                    if isKey(matchcomb_ids,key_tmp)
                        match_combinations{matchcomb_ids(key_tmp),1}(end+1) = bID1;
                    end        
                end
            end
            
            
            
            % tensordot must get int32
            legpos1 = int32(legpos1);
            legpos2 = int32(legpos2);
            freelegs1 = int32(freelegs1);
            freelegs2 = int32(freelegs2);
            
            for i = 1:(mID-1)
                if isempty(match_combinations{i,2}) || isempty(match_combinations{i,1})
                    continue;
                end
                for i1 = 1:length(match_combinations{i,1})
                    bID1 = match_combinations{i,1}(i1);
                    block1 = objvalues{bID1};
                    for i2 = 1:length(match_combinations{i,2})
                        bID2 = match_combinations{i,2}(i2);
                        block2 = othervalues{bID2};
                        naive_newkey = [objkeys{bID1},otherkeys{bID2}];
                        newkey = naive_newkey(newkey_layout);
                        if just_a_number
                            out = out + tensordot(block1,block2,legpos1,legpos2,freelegs1,freelegs2);
                        else
                            newblock = tensordot(block1,block2,legpos1,legpos2,freelegs1,freelegs2);
                            
                            if ~isKey(out.data, newkey)
                                out.data(newkey) = newblock;               
                            else
                                out.data(newkey) = out.data(newkey) + newblock;
                            end
                        end
                    end
                end
            end
end

