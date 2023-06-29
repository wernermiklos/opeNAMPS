function out = NTmatch_irreps(obj,irrep1,irrep2)
            % Matches two irreps, makes one irrep label of two. (Works only 
            % if the two irreps have equal values in every block, i.e. if the tensor is block diagonal)
            % irrep1 / irrep2 format: {'legname', depID}
            irrepIDs = NTlocate_irreps(obj, {irrep1,irrep2});
            if irrepIDs(1) < irrepIDs(2)
                irrep_kept = irrepIDs(1);
                irrep_disc = irrepIDs(2);
            elseif irrepIDs(2) < irrepIDs(1)
                irrep_kept = irrepIDs(2);
                irrep_disc = irrepIDs(1);
            else
                error('Error: the two irreps are the same: we cannot match them.')
            end
            newdeps = obj.dependencies;
            for legID = 1:length(newdeps)
                for depID = 1:length(newdeps{legID})
                    if newdeps{legID}(depID) == irrep_disc
                        newdeps{legID}(depID) = irrep_kept;
                    elseif newdeps{legID}(depID) > irrep_disc
                        newdeps{legID}(depID) = newdeps{legID}(depID) - 1;
                    end
                end
            end
            out = NAtensor(obj.leg_names,obj.leg_types,newdeps,obj.no_of_symmetries);
            for key = obj.data.keys
                if ~strcmp(key{1}((obj.no_of_symmetries*(irrep_kept-1)+1):(obj.no_of_symmetries*irrep_kept)), ...
                        key{1}((obj.no_of_symmetries*(irrep_disc-1)+1):(obj.no_of_symmetries*irrep_disc)))
                    error('Error: the two irreps have different values in some blocks')
                end
                newkey = key{1}([1:(obj.no_of_symmetries*(irrep_disc-1)),(obj.no_of_symmetries*irrep_disc + 1):end]);
                out.data(newkey) = obj.data(key{1});
            end
end

