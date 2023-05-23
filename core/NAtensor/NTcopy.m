function out = NTcopy( obj)
% Creates a (deep) copy of the object. (necessary because of the
%           containers.Map()) Using MatLab's new dictionary would make this
%           unnecessary.
            out = obj;   % this is just a shallow copy, map containers are not copied!!
            if obj.data.Count > 0
                out.data = containers.Map(obj.data.keys,obj.data.values,'UniformValues',false);
            else
                out.data = containers.Map('UniformValues',false);
            end
end

