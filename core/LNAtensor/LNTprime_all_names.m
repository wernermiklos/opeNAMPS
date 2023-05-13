function out = LNTprime_all_names(obj)
            % Puts a '~' at the and of all leg names.
            out = LNTcopy(obj);
            for legID = 1:length(obj.leg_names)
                out.leg_names{legID} = [out.leg_names{legID},'~'];
            end
            for symID = 1:obj.no_of_symmetries
                out.sub_tensors{symID} = NTprime_all_names(out.sub_tensors{symID});
            end


end

