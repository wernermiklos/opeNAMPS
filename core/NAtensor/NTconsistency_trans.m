function [other_leg_positions, other_irrep_positions, other_key_layout] = NTconsistency_trans(obj, other)
% Generates the consistency transformation from "other" to
% "obj" if the two tensors are consistent.
% ----------------------------------------
% other_leg_positions:     the leg positions in "other"
% other_irrep_positions:   the irrep positions in "other"
% key_layout:              the transformation rules for keys
%                          from other to obj
            if obj.no_of_symmetries ~= other.no_of_symmetries
                error('Error in NAtensor.consistency_trans(): the two NAtensors are not consistent.')
            end
            other_leg_positions = NTlocate_legs(other, obj.leg_names);
            if obj.irrep_number ~= other.irrep_number || length(obj.leg_names) ~= length(other.leg_names)
                error('Error in NAtensor.consistency_trans(): the two NAtensors are not consistent.')
            end
            other_irrep_positions = zeros(1,obj.irrep_number);
            for legID = 1:length(obj.leg_names)
                if ~strcmp(obj.leg_types{legID},other.leg_types{other_leg_positions(legID)})
                    error('Error in NAtensor.consistency_trans(): the two NAtensors are not consistent.')
                end
                if length(obj.dependencies{legID}) ~= length(other.dependencies{other_leg_positions(legID)})
                    error('Error in NAtensor.consistency_trans(): the two NAtensors are not consistent.')
                else    
                    for depID = 1:length(obj.dependencies{legID})
                        irrepID = obj.dependencies{legID}(depID);
                        other_irrepID = other.dependencies{other_leg_positions(legID)}(depID);
                        if other_irrep_positions(irrepID) == 0
                            other_irrep_positions(irrepID) = other_irrepID;
                        elseif other_irrep_positions(irrepID) ~= other_irrepID
                            error('Error in NAtensor.consistency_trans(): the two NAtensors are not consistent.')
                        end
                    end
                end
            end
            other_key_layout = zeros(1,obj.no_of_symmetries*obj.irrep_number);
            for irrepID = 1:obj.irrep_number
                other_key_layout((obj.no_of_symmetries*(irrepID-1)+1):(obj.no_of_symmetries*irrepID)) = ...
                   (obj.no_of_symmetries*(other_irrep_positions(irrepID)-1)+1):(obj.no_of_symmetries*other_irrep_positions(irrepID));
            end

end

