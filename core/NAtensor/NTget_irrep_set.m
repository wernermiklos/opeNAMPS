function out = NTget_irrep_set(obj,IrrepPosition)
            % Returns the set of possible (active) values of an irrep label specified by IrrepPosition.
            % ---
            % IrrepPosition:    specifies the irrep label. Format: {<legname>, depID}, where legname is the name of a
            %                   leg, while depID specifies the irrep label as the #th dependency of the leg.
            irrep_pos = NTlocate_irreps(obj,{IrrepPosition});
            key_layout = (obj.no_of_symmetries*(irrep_pos-1)+1):(obj.no_of_symmetries*irrep_pos);
            keylist = cell(1,obj.data.Count);
            i = 1;
            for key =obj.data.keys
                irrep_code = key{1}(key_layout);
                keylist{i} = irrep_code;
                i = i + 1;
            end
            keylist = unique(keylist);
            out = cell(1,length(keylist));
            for i = 1:length(out)
                out{i} = double(keylist{i});
            end
end

