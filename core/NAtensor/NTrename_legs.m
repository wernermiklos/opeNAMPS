function out = NTrename_legs(obj,oldnames,newnames)
    % Function to rename legs of an NAtensor.
    % ----
    % obj:        the NAtensor object
    % oldnames:   cell array of old names (can contain one or more leg names)
    % newnames:   cell array of new leg names (same length as oldnames!)
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

