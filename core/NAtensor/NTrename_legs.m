function out = NTrename_legs(obj,oldnames,newnames)
            legpos = NTlocate_legs(obj,oldnames);
            out = NTcopy(obj);
            for i = 1:length(legpos)
                out.leg_names{legpos(i)} = newnames{i};
            end
            if length(out.leg_names) ~= length(unique(out.leg_names))
                for i = 1:length(legpos)
                    out.leg_names{legpos(i)} = oldnames{i};
                end
                error('Error: new leg names are ambiguous.')
            end
end

