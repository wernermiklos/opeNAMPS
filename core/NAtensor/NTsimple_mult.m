function out = NTsimple_mult(obj, other, irrepmatch_rule)
            % Simple product of two tensors. (No legs are contracted)
            % ---
            % obj:     first NAtensor object
            % other:   second NAtensor object
            %       :::: Remark: all the leg names of the two NAtensors must be different.
            %       ::::         You may have to rename legs before calling this function.
            % 
            % irrepmatch_rule:  optional input. You may match one or more
            %                   irreps in the two tensors. 
            %                   format: {{{legname1_1, depID1_1},{legname2_1,
            % depID2_1}},{legname1_2, depID1_2},{legname2_2,depID2,2}}...}
            if nargin == 2
                irrepmatch_rule = {};
            elseif nargin == 3
            else
                error('NTsimple_mult has 2 or 3 arguments.')
            end
            if obj.no_of_symmetries ~= other.no_of_symmetries
                error('Error: the two tensors have different no_of_symmetries')
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
                legID1 = NTlocate_legs(obj,{irrepmatch_rule{i}{1}{1}});
                depID1 = irrepmatch_rule{i}{1}{2};
                legID2 = NTlocate_legs(other,{irrepmatch_rule{i}{2}{1}});
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
                irrepmatch_layout{i,1} = (obj.no_of_symmetries*(irrepID1-1)+1):(obj.no_of_symmetries*irrepID1);
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
            for key1 = obj.data.keys
                block1 = obj.data(key1{1});
                tmpblock1 = reshape(block1,[numel(block1),1]);
                tmpshape1 = [size(block1),ones(1,length(obj.leg_names)-ndims(block1))];
                for key2 = other.data.keys
                    kept = true;
                    for i = 1:length(irrepmatch_rule)
                        if ~strcmp(key1{1}(irrepmatch_layout{i,1}), key2{1}(irrepmatch_layout{i,2}))
                            kept = false;
                            break 
                        end
                    end
                    if kept
                        newkey_tmp = [key1{1},key2{1}];
                        newkey = newkey_tmp(newkey_layout);
                        block2 = other.data(key2{1});
                        tmpblock2 = reshape(block2,[1,numel(block2)]);
                        out.data(newkey) = reshape(tmpblock1*tmpblock2, ...
                            [tmpshape1,size(block2)]);
                    end
                end
            end


end

