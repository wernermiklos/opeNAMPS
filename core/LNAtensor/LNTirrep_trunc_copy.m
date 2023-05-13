function out = LNTirrep_trunc_copy(obj, KeptIrreps)
            % Returns a truncated copy of the tensor
            % KeptIrreps format: {{{'legname1',depID1},Irreplist1},...}
            % Irreps that are not specified are not truncated.
            % IMPORTANT!: This is the NAtensor compatibility version, i.e.
            % it is called in the same way as the NAtensor version.
            % REMARK: the result is not fully equivalent with the NAtensor
            % version. (see notes!)
            out = LNAtensor(obj.leg_names,obj.leg_types,obj.dependencies, obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                keptirreps_1sym = cell(1,length(KeptIrreps));
                for i = 1:length(KeptIrreps)
                    keptirreps_1sym{i} = {KeptIrreps{i}{1},zeros(1,length(KeptIrreps{i}{2}))};
                    for irrepID = 1:length(KeptIrreps{i}{2})
                        keptirreps_1sym{i}{2}(irrepID) = KeptIrreps{i}{2}{irrepID}(symID);
                    end
                    keptirreps_1sym{i}{2} = num2cell(unique(keptirreps_1sym{i}{2}));
                end
                out = LNTset_subtensor(out,symID,NTirrep_trunc_copy(obj.sub_tensors{symID},keptirreps_1sym));
            end


end

