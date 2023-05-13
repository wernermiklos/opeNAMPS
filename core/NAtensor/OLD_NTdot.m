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
            
            [irrep_classes,...
             irreps1_in_class,...
             irreps2_in_class,...
             just_a_number,...
             legpos1,...
             legpos2,...
             match_classes,...
             new_dependencies,...
             new_leg_names,...
             new_leg_types,...
             sum_up_classes] = OLD_NTdot_init(obj,other,legnames1,legnames2,renaming);
            if ~just_a_number
                out = NAtensor(new_leg_names, new_leg_types, new_dependencies, obj.no_of_symmetries);
            else
                out = 0.0;
            end
            
            % We create layouts to speed up things
            key1_matchpos = zeros(1,obj.no_of_symmetries*(length(sum_up_classes)+length(match_classes)));
            key2_matchpos = zeros(1,obj.no_of_symmetries*(length(sum_up_classes)+length(match_classes)));
            prefilter1 = false;    %If there is a match in tensor "obj", we will prefilter the blocks.
            prefilter2 = false;    %If there is a match in tensor "other", we will prefilter the blocks.
            filter1_testpos1 = [];     
            filter1_testpos2 = [];     
            filter2_testpos1 = [];     
            filter2_testpos2 = [];    
            j = 0;
            for classID = [sum_up_classes, match_classes]
                j = j + 1;
                irrep1_pos = irreps1_in_class{classID}(1);
                key1_matchpos((obj.no_of_symmetries*(j-1)+1):(obj.no_of_symmetries*(j))) = ...
                    (obj.no_of_symmetries*(irrep1_pos-1)+1):(obj.no_of_symmetries*(irrep1_pos));
                irrep2_pos = irreps2_in_class{classID}(1);
                key2_matchpos((other.no_of_symmetries*(j-1)+1):(other.no_of_symmetries*(j))) = ...
                    (other.no_of_symmetries*(irrep2_pos-1)+1):(other.no_of_symmetries*(irrep2_pos));
                if length(irreps1_in_class{classID}) > 1
                    prefilter1 = true;
                    for i = 1:(length(irreps1_in_class{classID})-1)
                        filter1_testpos1 = [filter1_testpos1, (obj.no_of_symmetries*((irreps1_in_class{classID}(i)-1)+1)):...
                                                                   (obj.no_of_symmetries*irreps1_in_class{classID}(i))];
                        filter1_testpos2 = [filter1_testpos2, (obj.no_of_symmetries*((irreps1_in_class{classID}(i+1)-1)+1)):...
                                                                   (obj.no_of_symmetries*irreps1_in_class{classID}(i+1))];
                    end
                end
                if length(irreps2_in_class{classID}) > 1
                    prefilter2 = true;
                    for i = 1:(length(irreps2_in_class{classID})-1)
                        filter2_testpos1 = [filter2_testpos1, (other.no_of_symmetries*((irreps2_in_class{classID}(i)-1)+1)):...
                                                                   (other.no_of_symmetries*irreps2_in_class{classID}(i))];
                        filter2_testpos2 = [filter2_testpos2, (other.no_of_symmetries*((irreps2_in_class{classID}(i+1)-1)+1)):...
                                                                   (other.no_of_symmetries*irreps2_in_class{classID}(i+1))];
                    end
                end
            end
                
            % Finally, before the contraction of blocks we sort the blocks. We create the dictionary MatchCombinations, in
            % which we collect block-key of self and other which blocks will be later contracted. The keys (key_tmp) of
            % MatchCombinations contain all the matched and sum_up irreps.
            matchcomb_ids = containers.Map();
            match_combinations = {};
            mID = 1;
            for key1 = obj.data.keys
                if prefilter1
                    if ~strcmp(key1{1}(filter1_testpos1),key1{1}(filter1_testpos2))
                        continue;
                    end
                end
                
                key_tmp = key1{1}(key1_matchpos);
                if isKey(matchcomb_ids,key_tmp)
                        match_combinations{matchcomb_ids(key_tmp)}{1}{end+1} = key1{1}; 
                else
                    matchcomb_ids(key_tmp) = mID;
                    mID = mID + 1;
                    match_combinations{end + 1} = {{key1{1}},{}};
                end      
            end
            
            for key2 = other.data.keys
                if prefilter2
                    if ~strcmp(key2{1}(filter2_testpos1),key2{1}(filter2_testpos2))
                        continue;
                    end
                end
                key_tmp = key2{1}(key2_matchpos);
                if isKey(matchcomb_ids,key_tmp)
                    match_combinations{matchcomb_ids(key_tmp)}{2}{end+1} = key2{1}; 
                end        
            end
            
           
            % We create a rule where the "new key" positions are in [key1,
            % key2]. This will be used to generate the new key from key1
            % and key2.
            newkey_positions = zeros(1,obj.no_of_symmetries*(length(irrep_classes)-length(sum_up_classes)));
            i = 1;
            for classID = 1:length(irrep_classes)
                if ~any(sum_up_classes == classID)
                    if irrep_classes{classID}{1}{1} == 1   % irrep class has member in obj
                        keypos = (obj.no_of_symmetries*(irrep_classes{classID}{1}{2} - 1) + 1):(obj.no_of_symmetries*(irrep_classes{classID}{1}{2}));
                        newkey_positions(i:(i+obj.no_of_symmetries-1)) = keypos;
                        i = i + obj.no_of_symmetries;
                    else
                        keypos = obj.no_of_symmetries*obj.irrep_number + ...
                            ((other.no_of_symmetries*(irrep_classes{classID}{1}{2} - 1) + 1):(other.no_of_symmetries*(irrep_classes{classID}{1}{2})));
                        newkey_positions(i:(i+other.no_of_symmetries-1)) = keypos;
                        i = i + obj.no_of_symmetries;
                    end
                end
            end
         
            
            % Finally we can perform block contractions.
            for key_tmp = matchcomb_ids.keys
                for key1 = match_combinations{matchcomb_ids(key_tmp{1})}{1}
                    shape1 = [size(obj.data(key1{1})),ones(1,length(obj.leg_names)-ndims(obj.data(key1{1})))];
                    for key2 = match_combinations{matchcomb_ids(key_tmp{1})}{2}
                        shape2 = [size(other.data(key2{1})),ones(1,length(other.leg_names)-ndims(other.data(key2{1})))];
                        totkey_tmp = [key1{1}, key2{1}];  % we concatenate the key1 and key2 strings.
                        if just_a_number
%                             out = out + NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
%                                           legpos1, legpos2,obj.shape(key1{1}), other.shape(key2{1}));
                              out = out + OLD_NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
                                           legpos1, legpos2,shape1, shape2);
                        else
                            newkey = totkey_tmp(newkey_positions);
                            if isKey(out.data,newkey)
%                                 out.data(newkey) = out.data(newkey) + ...
%                                                          NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
%                                                            legpos1, legpos2,obj.shape(key1{1}), other.shape(key2{1}));
                                out.data(newkey) = out.data(newkey) + ...
                                                         OLD_NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
                                                           legpos1, legpos2,shape1, shape2);
                            else
%                                 [out.data(newkey),out.shape(newkey)] = NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
%                                                             legpos1, legpos2,obj.shape(key1{1}), other.shape(key2{1}));
                                [out.data(newkey),~] = OLD_NTblockdot(obj.data(key1{1}),other.data(key2{1}), ...
                                                            legpos1, legpos2,shape1, shape2);
                            end
                        end
                    end
                end
            end      
end

