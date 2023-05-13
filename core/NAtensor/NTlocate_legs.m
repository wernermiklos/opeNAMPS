function out = NTlocate_legs(obj, LegNames)
% Returns the (internal) position of legs specified in LegNames
           out = zeros(1,length(LegNames));
           for i = 1:length(LegNames)
               pos = find(strcmp(obj.leg_names,LegNames{i}));
               if isempty(pos)
                  error(['Error in NAtensor.locate_legs(): leg "', LegNames{i}, '" is not found.']);
               end
               out(i) = pos;
           end 
end

