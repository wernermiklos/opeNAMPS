function out = NTmodify_irrep_values(obj,IrrepPosition,Gamma_to_Gammatild)
%NT_MODIFY_IRREP_VALUES Summary of this function goes here
%   Detailed explanation goes here
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

