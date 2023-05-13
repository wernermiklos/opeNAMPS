function out = NTcombine_legs(obj, LegNames, NewLegName)
%NTCOMBINE_LEGS Summary of this function goes here
%   Detailed explanation goes here
    legpos = NTlocate_legs(obj,LegNames);
    keptlegs = 1:length(obj.leg_names);
    keptlegs(legpos) = [];
    
    newlegnames = obj.leg_names;
    newlegnames(legpos) = [];
    newlegnames{end+1} = NewLegName;
    oldlegnum = length(obj.leg_names);

    comblegtype = obj.leg_types{legpos(1)};
    newlegtypes = obj.leg_types;
    newlegtypes(legpos) = [];
    newlegtypes{end+1} = comblegtype;

    newdeps = cell(1,length(newlegnames));
    for i = 1:(length(newlegnames)-1)
        newdeps{i} = obj.dependencies{keptlegs(i)};
    end

    for i = 1:length(legpos)
        if ~strcmp(obj.leg_types{legpos(i)},comblegtype)
            obj.leg_types{legpos}
            error('legtype mismatch in NTcombine_legs. All combined leg types must be the same!');
        end
        newdeps{end} = [newdeps{end}, obj.dependencies{legpos(i)}];
    end
    newdeps{end} = sort(unique(newdeps{end}));

    out = NAtensor(newlegnames,newlegtypes,newdeps,obj.no_of_symmetries);
    if oldlegnum > 1
        for key = obj.data.keys
            block = obj.data(key{1});
            tmpsize = [size(block),ones(1,oldlegnum-ndims(block))];
            newsize = [tmpsize(keptlegs),prod(tmpsize(legpos))];
            if ndims(newsize) == 1
                newsize = [newsize,1];
            end
            newblock = reshape(permute(obj.data(key{1}),[keptlegs,legpos]),newsize);
            out.data(key{1}) = newblock;
        end
    else
        error('Dont use NTcombine_legs for NAtensors with only one leg!!');
    end


end

