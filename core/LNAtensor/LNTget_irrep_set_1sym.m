function out = LNTget_irrep_set_1sym(obj,IrrepPosition)
            % Returns the set of possible (active) values of an irrep label specified by IrrepPosition.
            % ---
            % IrrepPosition:    specifies the irrep label. Format: {<legname>, depID}, where legname is the name of a
            %                   leg, while depID specifies the irrep label as the #th dependency of the leg.
            
            out = cell(1,obj.no_of_symmetries);
            for symID = 1:obj.no_of_symmetries
                out{symID} = NTget_irrep_set(obj.sub_tensors{symID},IrrepPosition);
            end


end

