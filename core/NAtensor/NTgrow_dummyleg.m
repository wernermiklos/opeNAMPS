function obj= NTgrow_dummyleg(obj, newlegname, newlegtype)
%NTGROW_LEG grows a dummy leg to the tensor obj. (sometimes maybe useful)
%   The dummy leg has no dependencies
% ---
% obj:            the NAtensor
% newlegname:     name of the new dummy leg
% newlegtype:     type of the new dummy leg.
obj = NTcopy(obj);
if any(cellfun(@(x) isequal(x,newlegname), obj.leg_names))
    error('redundantly defined leg names')
end
if newlegtype ~= 'i' && newlegtype ~= 'o'
    error('newlegtype should be "i" or "o".')
end
obj.leg_names = [obj.leg_names, newlegname];
obj.leg_types = [obj.leg_types, newlegtype];
obj.dependencies = [obj.dependencies, {[]}];

end

