function out = LNTget_irrep_set(obj,IrrepPosition)
            % Returns the set of possible (active) values of an irrep label specified by IrrepPosition.
            % ---
            % IrrepPosition:    specifies the irrep label. Format: {<legname>, depID}, where legname is the name of a
            %                   leg, while depID specifies the irrep label as the #th dependency of the leg.
            irrepsets_1sym = LNTget_irrep_set_1sym(obj,IrrepPosition);
            
            out = cell_cross_level1(irrepsets_1sym{1},irrepsets_1sym{2});
            for symID = 3:obj.no_of_symmetries
                out = cell_cross_level1(out,irrepsets_1sym{symID});
            end

end

