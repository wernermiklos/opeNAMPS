function out = NTget_block_shape(obj, IrrepPositions, IrrepIndices, LegNames)
            key = NTgenerate_key(obj,IrrepPositions,IrrepIndices);
            leg_positions = NTlocate_legs(obj,LegNames);
            if length(unique(leg_positions)) ~= length(obj.leg_names) || ...
                    length(leg_positions) ~= length(obj.leg_names)
                error(['Error in NAtensor.set_block(): LegNames does not ' ...
                    'specify unambigously every leg.'])
            end
            legnum  = length(obj.leg_names);
            if obj.data.isKey(key)
                tmpshape = [size(obj.data(key)),ones(1,legnum-ndims(obj.data(key)))];
                if legnum == 1
                    out = tmpshape(1);
                else
                    out = tmpshape(leg_positions);
                end
            else
                out = [];
            end

end

