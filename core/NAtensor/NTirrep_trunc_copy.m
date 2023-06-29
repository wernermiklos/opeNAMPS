function out = NTirrep_trunc_copy(obj, KeptIrreps)
            % Returns a truncated copy of the tensor. In the cell array
            % KeptIrreps we specify which quantum number sectors do we want
            % to keep. Irreps that are not listed in KeptIrreps are not truncated.
            % ----
            % KeptIrreps:   format: {{{'legname1',depID1},IrrepList1},...}
            %                    * {'legname1', depID1} specifies the irrep
            %                    * IrrepList1: cell array that contains the
            %                      irrep indices that are kept.
            % 
            out = NAtensor(obj.leg_names,obj.leg_types,obj.dependencies, obj.no_of_symmetries);
            irrepIDs = zeros(1,length(KeptIrreps));
            kept_irrepcodes = cell(1,length(KeptIrreps));
            for i = 1:length(KeptIrreps)
                legID = NTlocate_legs(obj,{KeptIrreps{i}{1}{1}});
                depID = KeptIrreps{i}{1}{2};
                irrepIDs(i) = obj.dependencies{legID}(depID);
                kept_irrepcodes{i} = cellfun(@char,KeptIrreps{i}{2},'UniformOutput',false);
            end
            for key = obj.data.keys
               kept = true;
               for i = 1:length(KeptIrreps)
                   irrepcode = key{1}((obj.no_of_symmetries*(irrepIDs(i)-1)+1):(obj.no_of_symmetries*irrepIDs(i)));
                   if ~any(strcmp(kept_irrepcodes{i},irrepcode))
                       kept = false;
                   end
               end
               if kept
                   out.data(key{1}) = obj.data(key{1});
               end
            end

end

