function out = NTneg(obj )
%Returns -obj
        out = struct();
        out.type = 'NAtensor';
        out.leg_names = obj.leg_names;
        out.leg_types = obj.leg_types;
        if ~isempty(obj.data)
            out.data = containers.Map(obj.data.keys,cellfun(@(x) -1*x, obj.data.values,'UniformOutput',false),'UniformValues',false);
        else
             out.data = containers.Map('UniformValues',false);
        end
        out.no_of_symmetries = obj.no_of_symmetries;
        out.irrep_number = obj.irrep_number;
        out.dependencies = obj.dependencies;
end

