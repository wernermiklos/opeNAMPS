function out = NTset_block(obj, IrrepPositions, IrrepIndices, LegNames, datablock)
            % Sets datablock as the block of the NAtensor object that is 
            % specified by IrrepPositions and IrrepIndices. LegNames
            % specifies the leg order of datablock.
            out = NTcopy(obj);
            key = NTgenerate_key(out,IrrepPositions,IrrepIndices);
            leg_positions = NTlocate_legs(out,LegNames);
            if length(unique(leg_positions)) ~= length(out.leg_names) || ...
                    length(leg_positions) ~= length(out.leg_names)
                error(['Error in NAtensor.set_block(): LegNames does not ' ...
                    'specify unambigously every leg.'])
            end
            if length(LegNames) > 1
                out.data(key) = ipermute(datablock,leg_positions);
            else
                out.data(key) = datablock(:);
            end
end
