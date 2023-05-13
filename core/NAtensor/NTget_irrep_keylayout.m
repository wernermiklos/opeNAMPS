function out = NTget_irrep_keylayout(obj,IrrepPosition)
%NTGET_IRREP_KEYLAYOUT Returns the key layout that corresponds to the irrep
%                      specified in IrrepPosition
% ---
% IrrepPosition:    specifies the irrep label. Format: {<legname>, depID}, where legname is the name of a
%                   leg, while depID specifies the irrep label as the #th dependency of the leg.
    irrep_pos = NTlocate_irreps(obj,{IrrepPosition});
    out = (obj.no_of_symmetries*(irrep_pos-1)+1):(obj.no_of_symmetries*irrep_pos);
end

