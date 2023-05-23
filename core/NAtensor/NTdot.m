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
            if verLessThan('matlab','9.12')
                error('MatLab R2022a or newer is required.')
            end

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
                
         
            dotinfo = NTdot_init(obj, other, legnames1, legnames2, renaming);

            matchkey_layout1 = zeros(1,obj.no_of_symmetries*length(dotinfo.matchirrep_list1));
            matchkey_layout2 = zeros(1,obj.no_of_symmetries*length(dotinfo.matchirrep_list2));
            for i = 1:length(dotinfo.matchirrep_list1)
                matchkey_layout1(obj.no_of_symmetries*(i-1)+1:obj.no_of_symmetries*i) = ...
                    obj.no_of_symmetries*(dotinfo.matchirrep_list1(i)-1)+1:obj.no_of_symmetries*dotinfo.matchirrep_list1(i);
                matchkey_layout2(other.no_of_symmetries*(i-1)+1:other.no_of_symmetries*i) = ...
                    other.no_of_symmetries*(dotinfo.matchirrep_list2(i)-1)+1:other.no_of_symmetries*dotinfo.matchirrep_list2(i);
            end
            
            
            % We determine how the new block key is generated from key1 and
            % key2
            outkey_layout1 = zeros(1,obj.irrep_number*obj.no_of_symmetries);
            outkey_layout2 = zeros(1,other.irrep_number*other.no_of_symmetries);
            for newirrepID = 1:dotinfo.new_irrep_number
                pos1 = find(dotinfo.new_irrepindices1 == newirrepID,1);
                if ~isempty(pos1)
                    outkey_layout1((obj.no_of_symmetries*(pos1-1)+1):(obj.no_of_symmetries*pos1)) = ...
                        (obj.no_of_symmetries*(newirrepID-1)+1):(obj.no_of_symmetries*newirrepID);
                else
                    pos2 = find(dotinfo.new_irrepindices2 == newirrepID,1);
                    outkey_layout2((other.no_of_symmetries*(pos2-1)+1):(other.no_of_symmetries*pos2)) = ...
                        (other.no_of_symmetries*(newirrepID-1)+1):(other.no_of_symmetries*newirrepID);
                end
            end
            
            outkey_length = obj.no_of_symmetries * dotinfo.new_irrep_number;
            objkeys = obj.data.keys;
            objvalues = obj.data.values;
            otherkeys = other.data.keys;
            othervalues = other.data.values;
            [keylist, tasklist] = NTdot_createtasklist(objkeys, ...
                                                       otherkeys, ...
                                                       int32(matchkey_layout1), ...
                                                       int32(matchkey_layout2), ...
                                                       int32(outkey_layout1), ...
                                                       int32(outkey_layout2), ...
                                                       int32(outkey_length));
            
            
            % we perform tasks
            legpos1 = dotinfo.legpos1;
            legpos2 = dotinfo.legpos2;                                     
            
            blocklist = cell(1,length(keylist));
            NumLegsObj = length(obj.leg_names);

            for taskID = 1:size(tasklist,1)
                if isempty(blocklist{tasklist(taskID,3)})
                    blocklist{tasklist(taskID,3)} = tensorprod(objvalues{tasklist(taskID,1)},...
                                                               othervalues{tasklist(taskID,2)}, ...
                                                               legpos1,...
                                                               legpos2,...
                                                               NumDimensionsA = NumLegsObj);
                else
                    blocklist{tasklist(taskID,3)} =  blocklist{tasklist(taskID,3)} + ...
                                                     tensorprod(objvalues{tasklist(taskID,1)},...
                                                                othervalues{tasklist(taskID,2)}, ...
                                                                legpos1,...
                                                                legpos2,...
                                                                NumDimensionsA = NumLegsObj);
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

