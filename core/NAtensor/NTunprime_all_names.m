function out =  NTunprime_all_names(obj)
%  Inverse of NTprime_all_names: removes one '~' from the end of all leg
%  names. Raises an error if not all names are primed originally.
% ---
% obj:    the NAtensor object
% ---
% out:    the NAtensor object, where the names are unprimed.
            out = NTcopy(obj);
            for legID = 1:length(out.leg_names)
                if out.leg_names{legID}(end) ~= '~'
                    error('Error: cannot unprime all the names. One or more of the names are not primed.')
                end
                out.leg_names{legID} = out.leg_names{legID}(1:(end-1));
            end

end

