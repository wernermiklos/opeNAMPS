function out = NTgrow_dummyleg(obj, newlegname, newlegtype)
%NTGROW_LEG Summary of this function goes here
%   Detailed explanation goes here
out = NTcopy(obj);
if any(cellfun(@(x) isequal(x,newlegname), obj.leg_names))
    error('redundantly defined leg names')
end
if newlegtype ~= 'i' && newlegtype ~= 'o'
    error('newlegtype should be "i" or "o".')
end
out.leg_names = [obj.leg_names, newlegname];
out.leg_types = [obj.leg_types, newlegtype];
out.dependencies = [obj.dependencies, {[]}];

end

