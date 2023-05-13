function out = LNTirrep_trunc_copy_LNA(obj, KeptIrreps)
            % Returns a truncated copy of the tensor
            % KeptIrreps format: {{{'legname1',depID1},Irreplist1},...}
            % Irreps that are not specified are not truncated.
            % IMPORTANT!: Irreplists are different from the NAtensor case.
            % - Irreplists are cell arrays of lenght no_of_symmetries.
            % - Each cell is the irreplist for the corresponding symmetry.
            % - If the cell is empty for a symmetry, that subtensor is not 
            %   truncated.  
            out = LNAtensor(obj.leg_names,obj.leg_types,obj.dependencies, obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                keptirreps_1sym = {};
                for i = 1:length(KeptIrreps)
                    if ~isempty(KeptIrreps{i}{2}{symID})
                        keptirreps_1sym{end+1} = {KeptIrreps{i}{1},KeptIrreps{i}{2}{symID}};
                    end
                end
                out = LNTset_subtensor(out,symID,NTirrep_trunc_copy(obj.sub_tensors{symID},keptirreps_1sym));
            end

end

