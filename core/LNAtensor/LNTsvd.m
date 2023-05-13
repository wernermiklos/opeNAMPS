function [left_tensor,schmidt_tensor,right_tensor] = LNTsvd(obj,left_leg_names,right_leg_names,bond_leg_names, bond_leg_types)
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
            
            left_sub_tensors = cell(1,obj.no_of_symmetries);
            right_sub_tensors = cell(1,obj.no_of_symmetries);
            schmidt_sub_tensors = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                [left_sub_tensors{symID},schmidt_sub_tensors{symID},right_sub_tensors{symID}] = NTsvd(obj.sub_tensors{symID},left_leg_names,right_leg_names,bond_leg_names,bond_leg_types);
            end
            left_tensor = LNAtensor(left_sub_tensors{1}.leg_names, ...
                                    left_sub_tensors{1}.leg_types, ...
                                    left_sub_tensors{1}.dependencies, ...
                                    obj.no_of_symmetries);
            right_tensor = LNAtensor(right_sub_tensors{1}.leg_names, ...
                                     right_sub_tensors{1}.leg_types, ...
                                     right_sub_tensors{1}.dependencies, ...
                                     obj.no_of_symmetries);
            schmidt_tensor = LNAtensor(schmidt_sub_tensors{1}.leg_names, ...
                                       schmidt_sub_tensors{1}.leg_types, ...
                                       schmidt_sub_tensors{1}.dependencies, ...
                                       obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
               left_tensor = LNTset_subtensor(left_tensor,symID,left_sub_tensors{symID});
               right_tensor = LNTset_subtensor(right_tensor,symID,right_sub_tensors{symID});
               schmidt_tensor = LNTset_subtensor(schmidt_tensor,symID,schmidt_sub_tensors{symID});
            end 


end

