function out = NTconj(obj)
% Generates the complex conjugate of the NAtensor. 
            % Remark: conjugation reverses the leg types, i.e. incoming
            % legs become outgoing ones, and outgoing legs become incoming
            % ones.
            if strcmp(obj.type,'LNAtensor')
                out = LNTconj(obj);
                return
            end
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
            
            out = struct();
            out.type = 'NAtensor';
            out.leg_names = obj.leg_names;
            out.leg_types = newlegtypes;
            if ~isempty(obj.data)
                out.data = containers.Map(obj.data.keys, cellfun(@conj,obj.data.values,'UniformOutput',false), 'UniformValues',false);
            else
                 out.data = containers.Map('UniformValues',false);
            end
            out.no_of_symmetries = obj.no_of_symmetries;
            out.irrep_number = obj.irrep_number;
            out.dependencies = obj.dependencies;
end

