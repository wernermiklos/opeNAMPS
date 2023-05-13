function out = LNTdiv(obj, val)
            % Overloaded function for obj  / val.
            if ~(isnumeric(val) && isequal(size(val),[1,1]))
                error('Error in LNAtensor / val : "val" must be 1D numerical variable.')
            end
            out = LNTmult(obj,(1.0 / val));
end


