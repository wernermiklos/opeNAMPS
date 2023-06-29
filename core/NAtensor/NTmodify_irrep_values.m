function out = NTmodify_irrep_values(obj,IrrepPosition,Gamma_to_Gammatild)
%NT_MODIFY_IRREP_VALUES Modifies the irrep values at IrrepPosition (one
%irrep can only be specified)
%---
% obj:                 the NAtensor
% IrrepPosition:       specifies the irrep whose indices we want to modify
%                      (one irrep can only be specified.)
%                      format:  {'legname', depID}
% Gamma_to_Gammatild:  array or function.
%                      newirrepindex = Gamma_to_Gammatild(oldirrepindex)
%                         * common example  Gamma_to_Gammatild can create the
%                                           conjugate representation indices
    out = obj;
    key_layout = NTget_irrep_keylayout(obj,IrrepPosition);
    datakeys = keys(obj.data);
    newkeys = cell(1,length(datakeys));
    datablocks = values(obj.data);

    for bID = 1:length(datakeys)   %If too slow then we should use a cellfun blackmagic.
        oldkey = datakeys{bID};
        newkey = oldkey;
        newkey(key_layout) = char(Gamma_to_Gammatild(double(oldkey(key_layout))));
        newkeys{bID} = newkey;
    end
    out.data = containers.Map(newkeys,datablocks,"UniformValues",false);

end

