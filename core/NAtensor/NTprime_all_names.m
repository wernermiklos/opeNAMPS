function out = NTprime_all_names(obj)
            % Puts a '~' at the and of all leg names.
            out = NTcopy(obj);
            for legID = 1:length(obj.leg_names)
                out.leg_names{legID} = [out.leg_names{legID},'~'];
            end
end

