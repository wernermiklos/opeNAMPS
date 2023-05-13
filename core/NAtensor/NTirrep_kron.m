function out = NTirrep_kron(obj, other, varargin)
            % Returns irrep_kron of obj and other.
            % If as third argument the string 'LNAtensor' is passed, the
            % output is LNAtensor, otherwise the full NAtensor is produced.
            if ~strcmp(other.type,'NAtensor')
                error('Error: other must be an NAtensor too.')
            end
            legnum = length(obj.leg_names);
            if ~isempty(varargin)
                if strcmp(varargin{1},'LNAtensor')
                    % warning('LNAtensors are temporarily not used.');
                     if (obj.no_of_symmetries ~= 1) || (other.no_of_symmetries ~= 1)
                         warning('LNAtensor output is possible only if no_of_symmetries is 1 for each tensor')
                     else
                         out = LNAtensor(obj.leg_names, obj.leg_types, obj.dependencies, 2);
                         out = LNTset_subtensor(out, 1, obj);
                         out = LNTset_subtensor(out, 2, other);
                         return
                     end
                else
                    warning([varargin{1}, ' argument is not accepted. Try LNAtensor instead, if you want an LNAtensor output.'])
                end
            end
            other_leg_positions = NTlocate_legs(other,obj.leg_names);
            if obj.irrep_number ~= other.irrep_number || length(obj.leg_names) ~= length(other.leg_names)
                error('Error: the two NAtensors are not consistent.')
            end
            other_irrep_positions = zeros(1,obj.irrep_number);
            for legID = 1:length(obj.leg_names)
                if ~strcmp(obj.leg_types{legID},other.leg_types{other_leg_positions(legID)})
                    error('Error: the two NAtensors are not consistent.')
                end
                if length(obj.dependencies{legID}) ~= length(other.dependencies{other_leg_positions(legID)})
                    error('Error: the two NAtensors are not consistent.')
                else    
                    for depID = 1:length(obj.dependencies{legID})
                        irrepID = obj.dependencies{legID}(depID);
                        other_irrepID = other.dependencies{other_leg_positions(legID)}(depID);
                        if other_irrep_positions(irrepID) == 0
                            other_irrep_positions(irrepID) = other_irrepID;
                        elseif other_irrep_positions(irrepID) ~= other_irrepID
                            error('Error: the two NAtensors are not consistent.')
                        end
                    end
                end
            end
            other_key_layout = zeros(1,other.no_of_symmetries*other.irrep_number);
            for irrepID = 1:other.irrep_number
                other_key_layout((other.no_of_symmetries*(irrepID-1)+1):(other.no_of_symmetries*irrepID)) = ...
                   (other.no_of_symmetries*(other_irrep_positions(irrepID)-1)+1):(other.no_of_symmetries*other_irrep_positions(irrepID));
            end
            out = NAtensor(obj.leg_names,obj.leg_types,obj.dependencies,obj.no_of_symmetries + other.no_of_symmetries);
            for key1cell = obj.data.keys
                key1 = key1cell{1};
                block1 = obj.data(key1);
                tmpshape1 = [size(block1),ones(1,legnum-ndims(block1))];
                for key2cell = other.data.keys
                    key2 = key2cell{1};
                    block2 = other.data(key2);
                    newkey = '';
                    for irrepID = 1:obj.irrep_number
                        newkey = [newkey, key1((obj.no_of_symmetries*(irrepID-1)+1):(obj.no_of_symmetries*irrepID))];
                        newkey = [newkey, key2(other_key_layout((other.no_of_symmetries*(irrepID-1)+1):(other.no_of_symmetries*irrepID)))];
                    end
                    tmpshape2 = [size(block2),ones(1,legnum-ndims(block2))];
                    [out.data(newkey), ~] = NTblock_kron(block1,...
                                                                         permute(block2, other_leg_positions),...
                                                                         tmpshape1,...
                                                                         tmpshape2(other_leg_positions));
                end
            end

end

