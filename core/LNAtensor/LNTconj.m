function out = LNTconj(obj)
            % Generates the complex conjugate of the LNAtensor. 
            % Remark: conjugation reverses the leg types, i.e. incoming
            % legs become outgoing ones, and outgoing legs become incoming
            % ones.
            newlegtypes = cell(1,length(obj.leg_types));
            for i = 1:length(obj.leg_types)
                if obj.leg_types{i} == 'i'
                    newlegtypes{i}='o';
                elseif obj.leg_types{i} == 'o'
                    newlegtypes{i}='i';
                else
                    error('Error in NAtensor.leg_types.')
                end
            end
            out = LNAtensor(obj.leg_names, newlegtypes, obj.dependencies, obj.no_of_symmetries);
            out.sub_just_ones = obj.sub_just_ones;
            out.sub_just_nums = obj.sub_just_nums;
            for symID = 1:obj.no_of_symmetries
                out.sub_tensors{symID} = NTconj(obj.sub_tensors{symID});
            end

end

