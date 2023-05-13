function out = NTdiv(obj,val)
%NTDIV Returns obj / val.
     if ~(isnumeric(val) && isequal(size(val),[1,1]))
                error('Error in LNAtensor / val : "val" must be 1D numerical variable.')
     end
     out = NTmult(obj,1.0 / val);

end

