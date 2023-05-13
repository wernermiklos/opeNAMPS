function out = LNTrename_legs(obj,oldnames,newnames)
            % Renames legs that are listed in oldnames. 
            % The new names are listed in newnames.
            % Legnames that are not listed in oldnames are kept.
            out = LNTcopy(obj);
            for symID = 1:obj.no_of_symmetries
                out.sub_tensors{symID} = NTrename_legs(obj.sub_tensors{symID},oldnames,newnames);
            end
            out.leg_names = out.sub_tensors{1}.leg_names;  % We change also the objects leg_names list.


end

