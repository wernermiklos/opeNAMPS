function out = NTmult(obj,val)
%NTMULT returns obj*val
    out = struct();
    out.type = 'NAtensor';
    out.leg_names = obj.leg_names;
    out.leg_types = obj.leg_types;
    tmpkeys = obj.data.keys;
    tmpvalues = obj.data.values;
    tmpvalues = cellfun(@(x) val*x, tmpvalues , 'UniformOutput', false);
    if ~isempty(tmpkeys)
        out.data = containers.Map(tmpkeys,tmpvalues,'UniformValues',false);
    else
        out.data = containers.Map('UniformValues',false);
    end
    out.no_of_symmetries = obj.no_of_symmetries;
    out.irrep_number = obj.irrep_number;
    out.dependencies = obj.dependencies;
end

