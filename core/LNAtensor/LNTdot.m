function out = LNTdot(a, b, legnames_a, legnames_b, renaming_ab)
            % Calculates dot product (leg contraction) of "self" and "other". The contracted legs are specified in
            % "LegNames1" and "LegNames2". In Renaming1 and Renaming2 one can specify new names of (remaining) legs after
            % the contraction. If a leg is not listed in Renaming1 or Renaming2, the old name will be used.
            % If the list LegNames1 and LegNames2 is empty, NAtensor.Simplemult() will be called (Renamings are performed!)
            % ----
            % a:            NAtensor or LNAtensor
            % b:            NAtensor or LNAtensor
            % legnames_a         list of leg names of "a" that are contracted.
            % legnames_b         list of leg names of "b" that are contracted with the legs of "obj". Length of
            %                   legnames_a must equal length of legnames_b. Contractions are performed pairwise (i.e.
            %                   legnames_a{0} is contracted with legnames_b{0}, ...). The contracted legs must be of the same
            %                   kind, (same dependency structure, same dimensions in the two tensors).
            % renaming          List of leg renamings in "a" and "b". 
            %                   renaming{1} contains the renamings of "obj"
            %                   renaming{2} containes the renamings of
            %                   "other"
            %                   Format: {{{'oldname1_1', 'newname1_1'}, ...
            %                             {'oldname1_2', 'newname1_2'}, ...},
            %                            {{'oldname2_1', 'newname2_1'}, ...
            %                             {'oldname2_2', 'newname2_2'}}}
            
            if verLessThan('matlab','9.12')
                actual_tensordot = @tensordot_C;   % for release < R2022a we must use the C-MEX version
            else
                actual_tensordot = @tensordot;     % for release <= R2022a we use the wrapped built in version.
            end
            
            if nargin == 4
                renaming_ab = {{},{}};
            elseif nargin == 5
            else
                error('LNTdot() has 4 or 5 parameters. The 4th (optional) parameter is "renaming"')
            end
            if a.no_of_symmetries ~= b.no_of_symmetries
                    error('No_of_symmetries must be the same for the contracted tensors!')
            end
            if strcmp(a.type,'LNAtensor') && strcmp(b.type,'LNAtensor')   % If both tensors are LNAtensors, we have an easy job:
                new_subtensors = cell(1,a.no_of_symmetries);
                for symID = 1:a.no_of_symmetries
                    new_subtensors{symID} = NTdot(a.sub_tensors{symID},b.sub_tensors{symID},legnames_a,legnames_b,renaming_ab);
                end
                if isnumeric(new_subtensors{1})  % If the result is just a number (with no legs)
                     out = prod(cell2mat(new_subtensors));
                else
                    out = LNAtensor(new_subtensors{1}.leg_names, new_subtensors{1}.leg_types, new_subtensors{1}.dependencies, a.no_of_symmetries);
                    for symID = 1:a.no_of_symmetries
                        out = LNTset_subtensor(out,symID,new_subtensors{symID});
                    end
                end
                return
            elseif strcmp(a.type,'LNAtensor') && strcmp(b.type,'NAtensor')
                obj = a;
                other = b;
                legnames1 = legnames_a;
                legnames2 = legnames_b;
                renaming = renaming_ab;
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'LNAtensor')
                obj = b;
                other = a;
                legnames1 = legnames_b;
                legnames2 = legnames_a;
                renaming = {renaming_ab{2},renaming_ab{1}};     % We switch the two cells.
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'NAtensor')
                out = NTdot(a,b,legnames_a,legnames_b,renaming_ab);
                return
            else
                error('Woops');
            end
            
            % othersub_dummy = NAtensor(other.leg_names, other.leg_types, other.dependencies, 1);
            
            dotinfo = NTdot_init(obj,other,legnames1,legnames2,renaming);
        
            matchkey_layout_list1 = cellfun(@(x) zeros(1,length(dotinfo.matchirrep_list1)),cell(1,obj.no_of_symmetries),'UniformOutput',false);
            matchkey_layout_list2 = cellfun(@(x) zeros(1,length(dotinfo.matchirrep_list2)),cell(1,other.no_of_symmetries), 'UniformOutput',false);
            for symID = 1:obj.no_of_symmetries
                matchkey_layout_list1{symID} = dotinfo.matchirrep_list1;
                matchkey_layout_list2{symID} = other.no_of_symmetries*(dotinfo.matchirrep_list2-1)+symID;
            end
            
            outkey_layout_list1 = cellfun(@(x) zeros(1,obj.irrep_number),cell(1,obj.no_of_symmetries),'UniformOutput',false);
            outkey_layout2 = zeros(1,other.irrep_number*other.no_of_symmetries);
            for newirrepID = 1:dotinfo.new_irrep_number
                pos1 = find(dotinfo.new_irrepindices1 == newirrepID,1);
                if ~isempty(pos1)
                    for symID = 1:obj.no_of_symmetries
                        outkey_layout_list1{symID}(pos1) = ...
                        obj.no_of_symmetries*(newirrepID-1)+symID;
                    end
                else
                    pos2 = find(dotinfo.new_irrepindices2 == newirrepID,1);
                    outkey_layout2((other.no_of_symmetries*(pos2-1)+1):(other.no_of_symmetries*pos2)) = ...
                        (other.no_of_symmetries*(newirrepID-1)+1):(other.no_of_symmetries*newirrepID);
                end
            end
            outkey_length = obj.no_of_symmetries * dotinfo.new_irrep_number;
            
            
            if ~dotinfo.just_a_number
                out = NAtensor(dotinfo.new_leg_names, dotinfo.new_leg_types, dotinfo.new_dependencies, obj.no_of_symmetries);
            else
                out = 0.0;
            end
            
            objkey_lists = cellfun(@(x) keys(x.data), obj.sub_tensors,'UniformOutput',false);
            otherkeys = keys(other.data);
            
            
%             [keylist, tasklist] = LNTdot_createtasklist(objkey_lists, ...
%                                                         otherkeys,...
%                                                         matchkey_layout_list1,...
%                                                         matchkey_layout_list2,...
%                                                         outkey_layout_list1,...
%                                                         outkey_layout2,...
%                                                         outkey_length,...
%                                                         obj.no_of_symmetries);

            %We call the C-version
            matchkey_layout_list1_int32 = cell(1,obj.no_of_symmetries);
            matchkey_layout_list2_int32 = cell(1,obj.no_of_symmetries);
            outkey_layout_list1_int32 = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                matchkey_layout_list1_int32{symID} = int32(matchkey_layout_list1{symID});
                matchkey_layout_list2_int32{symID} = int32(matchkey_layout_list2{symID});
                outkey_layout_list1_int32{symID} = int32(outkey_layout_list1{symID});
            end
            [keylist, tasklist] = LNTdotNT_createtasklist(objkey_lists, ...                                             
                                                        otherkeys,...
                                                        matchkey_layout_list1_int32,...
                                                        matchkey_layout_list2_int32,...
                                                        outkey_layout_list1_int32,...
                                                        int32(outkey_layout2),...
                                                        int32(outkey_length),...
                                                        int32(obj.no_of_symmetries));
            
            
            % Finally, before the contraction of blocks we sort the blocks. We create the dictionary MatchCombinations, in
            % which we collect block-key of self and other which blocks will be later contracted. The keys (key_tmp) of
            % MatchCombinations contain all the matched and sum_up irreps.
            
            legpos1 = int32(dotinfo.legpos1);
            legpos2 = int32(dotinfo.legpos2);
            keptlegs1 = int32(dotinfo.keptlegs1);
            keptlegs2 = int32(dotinfo.keptlegs2);                                       
            blocklist = cell(1,length(keylist));

            sub_blocks = cellfun(@(x) values(x.data), obj.sub_tensors,'UniformOutput',false);
            other_blocks = values(other.data);
            sub_bIDs = zeros(1,obj.no_of_symmetries);
            tmp_blocks = cell(1,obj.no_of_symmetries);
            outlegnum = length(dotinfo.new_leg_names);
            aux_outleg_list = 1:outlegnum;

            for taskID = 1:size(tasklist,1)
                for symID = 1:obj.no_of_symmetries
                    if tasklist(taskID,symID)==sub_bIDs(symID)
                        continue
                    else
                        for i = symID:obj.no_of_symmetries
                            sub_bIDs(i) = tasklist(taskID,i);
                            tmp_subblock = sub_blocks{i}{tasklist(taskID,i)};
                            if i == 1
                                tmp_blocks{i} = tmp_subblock;
                            else
                                if obj.sub_just_ones(i)
                                    tmp_blocks{i} = tmp_blocks{i-1};
                                elseif obj.sub_just_nums(i)
                                    tmp_blocks{i} = tmp_subblock*tmp_blocks{i-1};
                                else
                                    tmp_blocks{i} = NTblock_kron(tmp_blocks{i-1},tmp_subblock,...
                                                                 size(tmp_blocks{i-1},aux_outleg_list), ...
                                                                 size(tmp_subblock,aux_outleg_list));
                                    
                                end
                            end
                        end
                        break;
                    end
                end
                if isempty(blocklist{tasklist(taskID,end)})
                    blocklist{tasklist(taskID,end)} = actual_tensordot(tmp_blocks{end},...
                                                                other_blocks{tasklist(taskID,end-1)}, ...
                                                                legpos1,legpos2,keptlegs1,keptlegs2);
                else
                    blocklist{tasklist(taskID,end)} =  blocklist{tasklist(taskID,end)} + ...
                                                          actual_tensordot(tmp_blocks{end},...
                                                                    other_blocks{tasklist(taskID,end-1)}, ...
                                                                    legpos1,legpos2,keptlegs1,keptlegs2);
                end
            end
            
            % we create the out variable
            if ~dotinfo.just_a_number
                out = NAtensor(dotinfo.new_leg_names, dotinfo.new_leg_types, dotinfo.new_dependencies, obj.no_of_symmetries);
                if ~isempty(keylist)
                    out.data = containers.Map(keylist,blocklist,'UniformValues',false);
                end     
            else
                out = 0;
                if ~isempty(keylist)
                    out = blocklist{1};
                end
            end



end

