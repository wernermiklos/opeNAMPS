function out = NTsplit_irrep(obj,IrrepPos,legs_dep_on_new)
            % Splits one irrep label into two. (Usually used for splitting
            % the irrep label of irrep-'diagonal' tensors.
            %     :::: this is the inverse of NTmatch_irreps()
            % ---
            % IrrepPos:            {legname,depID}
            % legs_dep_on_new:     cell array leg names that depend on the
            %                      new irrep ID instead of the old one.

            old_irrepID = NTlocate_irreps(obj,{IrrepPos});
            legIDs = NTlocate_legs(obj,legs_dep_on_new);
            newdeps = obj.dependencies;
            for legID = legIDs 
                for depID = 1:length(obj.dependencies{legID})
                    if newdeps{legID}(depID) == old_irrepID
                        newdeps{legID}(depID) = obj.irrep_number + 1;
                    end
                end
            end
            
            out = NAtensor(obj.leg_names,obj.leg_types,newdeps,obj.no_of_symmetries);
            for key_cell = obj.data.keys
                old_key = key_cell{1};
                splitted_irrep = old_key(((old_irrepID-1)*obj.no_of_symmetries+1):(old_irrepID*obj.no_of_symmetries));
                new_key = [old_key,splitted_irrep];
                out.data(new_key) = obj.data(old_key);
            end


end

