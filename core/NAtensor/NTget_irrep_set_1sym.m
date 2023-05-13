function out = NTget_irrep_set_1sym(obj,IrrepPosition)
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
            keylist_1sym = cell(1,obj.no_of_symmetries);
            out = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                keylist_1sym{symID} = cell(1,length(keylist));
                for i = 1:length(keylist)
                    keylist_1sym{symID}{i} = keylist{i}(symID);
                end
                keylist_1sym{symID} = unique(keylist_1sym{symID});
                out{symID} = cell(1,length(keylist_1sym{symID}));
                for i = 1:length(out{symID})
                    out{symID}{i} = double(keylist_1sym{symID}{i});
                end
            end
end

