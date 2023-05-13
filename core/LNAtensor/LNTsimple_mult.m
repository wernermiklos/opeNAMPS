function out = LNTsimple_mult(a, b, irrepmatch_rule_ab)
            % Simple product of two tensors. (No legs are contracted)
            % irrepmatch_rule format: {{{legname1_1, depID1_1},{legname2_1,
            % depID2_1}},{legname1_2, depID1_2},{legname2_2,depID2,2}}...}
            if nargin == 2
                irrepmatch_rule_ab = {};
            elseif nargin == 3
            else
                error('simple_mult has 2 or 3 arguments.')
            end
            
            if a.no_of_symmetries ~= b.no_of_symmetries
                error('Error: the two tensors have different no_of_symmetries')
            end
            if strcmp(a.type,'LNAtensor') && strcmp(b.type,'LNAtensor')
                new_sub_tensors = cell(1,a.no_of_symmetries);
                for symID = 1:a.no_of_symmetries
                    new_sub_tensors{symID} = NTsimple_mult(a.sub_tensors{symID},b.sub_tensors{symID},irrepmatch_rule_ab);
                end
                out = LNAtensor(new_sub_tensors{1}.leg_names, ...
                                new_sub_tensors{1}.leg_types, ...
                                new_sub_tensors{1}.dependencies, ...
                                a.no_of_symmetries);
                for symID = 1:a.no_of_symmetries
                    out = LNTset_subtensor(out,symID,new_sub_tensors{symID});
                end
                return
            elseif strcmp(a.type,'LNAtensor') && strcmp(b.type,'NAtensor')
                obj = a;
                other = b;
                irrepmatch_rule = irrepmatch_rule_ab;
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'LNAtensor')
                obj = b;
                other = a;
                irrepmatch_rule = cell(1,length(irrepmatch_rule_ab));
                for i = 1:length(irrepmatch_rule)
                    irrepmatch_rule{i} = {irrepmatch_rule_ab{i}{2}, irrepmatch_rule_ab{i}{1}};
                end
            end
            
            newlegnames = [obj.leg_names, other.leg_names];
            if length(newlegnames) ~= length(unique(newlegnames))
                error('Error: amibguously defined new leg names.')
            end
            newlegtypes = [obj.leg_types, other.leg_types];
            
            % We create first the naiive newdependencies (irrepmatch_rule
            % is ignored first).
            newdependencies = cell(1,length(newlegnames));
            for legID = 1:length(obj.leg_names)
                newdependencies{legID} = obj.dependencies{legID};
            end
            for legID = (length(obj.leg_names)+1):(length(obj.leg_names) + length(other.leg_names))
                legID2 = legID - length(obj.leg_names);
                newdependencies{legID} = other.dependencies{legID2} + obj.irrep_number;  %We shift the irrepID's
            end
                  
            % Now we "decode" irrepmatch_rule, and create new irrepIDs
            newirrepIDs = 1:(obj.irrep_number + other.irrep_number);
            irrepmatch_layout = cell(length(irrepmatch_rule),2);
            for i = 1:length(irrepmatch_rule)
                legID1 = LNTlocate_legs(obj,{irrepmatch_rule{i}{1}{1}});
                depID1 = irrepmatch_rule{i}{1}{2};
                legID2 = LNTlocate_legs(other,{irrepmatch_rule{i}{2}{1}});
                depID2 = irrepmatch_rule{i}{2}{2};
                irrepID1 = obj.dependencies{legID1}(depID1);
                irrepID2 = other.dependencies{legID2}(depID2);
                tmpID1 = newirrepIDs(irrepID1);
                tmpID2 = newirrepIDs(irrepID2 + obj.irrep_number);
                smaller = min([tmpID1,tmpID2]);
                larger = max([tmpID1,tmpID2]);
                if smaller ~= larger
                    for j = 1:length(newirrepIDs)
                        if newirrepIDs(j) == larger
                            newirrepIDs(j) = smaller;
                        elseif newirrepIDs(j) > larger
                            newirrepIDs(j) = newirrepIDs(j) - 1;
                        end
                    end
                end
                irrepmatch_layout{i,1} = irrepID1;
                irrepmatch_layout{i,2} = (other.no_of_symmetries*(irrepID2-1)+1):(other.no_of_symmetries*irrepID2);
            end
            
            % We transform the newdependencies list according to
            % newirrepIDs
            for legID = 1:length(newdependencies)
                newdependencies{legID} = newirrepIDs(newdependencies{legID});
            end
            
            out = NAtensor(newlegnames, newlegtypes, newdependencies, obj.no_of_symmetries);
            
            % We determine how the new key will be created from key1 and
            % key2.
            newkey_layout = zeros(1,max(newirrepIDs)*obj.no_of_symmetries);
            pos = 0;
            for i = 1:length(newirrepIDs)
                if newirrepIDs(i) > pos
                    pos = pos + 1;
                    newkey_layout((obj.no_of_symmetries*(pos-1)+1):(obj.no_of_symmetries*pos)) = ...
                        (obj.no_of_symmetries*(i-1)+1):(obj.no_of_symmetries*i);
                end
            end
            
            % Finally we perform the multiplications
            
            for key2 = other.data.keys
                kept_subkeys = cell(1,obj.no_of_symmetries);
                for symID = 1:obj.no_of_symmetries
                    for subkey1 = obj.sub_tensors{symID}.data.keys
                        kept = true;
                        for i = 1:length(irrepmatch_rule)
                            if ~strcmp(subkey1{1}(irrepmatch_layout{i,1}), key2{1}(irrepmatch_layout{i,2}(symID)))
                                kept = false;
                                break 
                            end
                        end
                        if kept
                            kept_subkeys{symID}{end+1} = subkey1{1}; 
                        end
                    end
                end
                tmp = cell(1,length(kept_subkeys{1})*length(kept_subkeys{2}));
                i = 1;
                for subkey1 = kept_subkeys{1}
                    for subkey2 = kept_subkeys{2}
                        tmp{i} = {subkey1{1}, subkey2{1}};
                        i = i + 1;
                    end
                end
                for symID = 3:obj.no_of_symmetries
                    tmp_old = tmp;
                    tmp = cell(1,length(tmp_old)*length(kept_subkeys{symID}));
                    i = 1;
                    for j = 1:length(tmp_old)
                        for subkey = kept_subkeys{symID}
                            tmp{i} = [tmp_old{j},{subkey{1}}];
                            i = i + 1;
                        end
                    end    
                end
                subkey_sets = tmp;
                for subkeys = subkey_sets
                    newkey_tmp = '';
                    for irrepID = 1:obj.irrep_number
                        for symID = 1:obj.no_of_symmetries
                            newkey_tmp = [newkey_tmp, subkeys{1}{symID}(irrepID)];
                        end
                    end
                    newkey_tmp = [newkey_tmp, key2{1}];  % we concatenate the key1 and key2 strings.
                    newkey = newkey_tmp(newkey_layout);
                    block1 = [1];
                    shape1 = ones(1,length(obj.leg_names));
                    for symID = 1:obj.no_of_symmetries
                        if obj.sub_just_ones(symID)
                        elseif obj.sub_just_nums(symID)
                            subblock = obj.sub_tensors{symID}.data(subkeys{1}{symID});
                            block1 = block1 * subblock.values;
                        else
                            block_tmp =  obj.sub_tensors{symID}.data(subkeys{1}{symID});
                            shape_tmp = size(block_tmp);
                            shape_tmp = [shape_tmp, ones(1,length(other.leg_names)-length(shape_tmp))];
                            [block1,shape1] = NTblock_kron(block1,block_tmp,shape1,shape_tmp);
                        end
                    end
                    block2 = other.data(key2{1});
                    shape2 = size(block2);
                    out.data(newkey) = reshape(block1(:)*(block2(:).'), ...
                            [shape1,shape2]);
                end
            end

end

