function out = NTget_block(obj, IrrepPositions, IrrepIndices, LegNames)
            % Returns a single block of an NAtensor.
            % ---
            % - IrrepPositions: format {{legname1,depID1},{legname2,depID2},...}
            % - IrrepIndices: cell array of irrep quantum numbers whose
            %                 order is specified by IrrepPositions
            % - LegNames: list af (all) leg names that specifies how the
            %             indices of the output are permuted.
            key = NTgenerate_key(obj,IrrepPositions,IrrepIndices);
            leg_positions = NTlocate_legs(obj,LegNames);
            if length(unique(leg_positions)) ~= length(obj.leg_names) || ...
                    length(leg_positions) ~= length(obj.leg_names)
                error(['Error in NAtensor.set_block(): LegNames does not ' ...
                    'specify unambigously every leg.'])
            end
            if obj.data.isKey(key)
                out = permute(obj.data(key),leg_positions);
            else
                out = [];
            end
end

