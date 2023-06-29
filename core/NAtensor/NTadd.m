function out = NTadd(obj, other)
           % Adds two NAtensors ("obj" and "other"),  if their structure 
           % (leg_names, leg_types, dependency structure, block shapes) are consistent.
           %
           % It is also allowed to add a tensor to the 0 number, then the
           % nonzero tensor is just copied.
           % ----
           % obj:    - NAtensor
           % other:  - the other NAtensor
           if isequal(other,0)
               out = NTcopy(obj);
           elseif isequal(obj,0)
               out = NTcopy(other);
           else
                [other_leg_positions, ~, other_key_layout] = ...
                    NTconsistency_trans(obj, other);
                out = NTcopy(obj);
                for key = other.data.keys
                    newkey = key{1}(other_key_layout);
                    if out.data.isKey(newkey)
                        out.data(newkey) = out.data(newkey) + ...
                            permute(other.data(key{1}),other_leg_positions);
                    else
                        out.data(newkey) = permute(other.data(key{1}),other_leg_positions);
                    end
                end
                
            end
end

