function out =  NTunprime_all_names(obj)
            out = NTcopy(obj);
            for legID = 1:length(out.leg_names)
                if out.leg_names{legID}(end) ~= '~'
                    error('Error: cannot unprime all the names. One or more of the names are not primed.')
                end
                out.leg_names{legID} = out.leg_names{legID}(1:(end-1));
            end

end

