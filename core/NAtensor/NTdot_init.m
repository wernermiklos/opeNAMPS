function out = NTdot_init(obj, other, legnames1, legnames2, renaming)
%NTDOT_INIT init function called only by NTdot and LNTdot.
    if ~strcmp(other.type,'NAtensor') && ~strcmp(other.type,'LNAtensor')
        error('Error in NTdot(): other must be an NAtensor or LNAtensor.')
    end

    if length(legnames1) ~= length(legnames2)
        error('The number of multiplied legs are not the same in the two tensor!')
    end
    if obj.no_of_symmetries ~= other.no_of_symmetries
        error('The number of symmetries of the two tensors are not equal')
    end
    legpos1 = NTlocate_legs(obj, legnames1);
    legpos2 = NTlocate_legs(other, legnames2);
    
    if length(unique(legpos1)) ~= length(legpos1) || ...
           length(unique(legpos2)) ~= length(legpos2)
        error('Repetition in legnames1 or legnames2')
    end

    objlegtypes_tmp = obj.leg_types(legpos1);
    otherlegtypes_tmp = other.leg_types(legpos2);
    if any(cellfun(@(x,y) x==y, objlegtypes_tmp, otherlegtypes_tmp))
        error('Legs with different types ("i","o") can only be contracted.')
    end
    
    keptlegs1 = 1:length(obj.leg_names);
    keptlegs1(legpos1) = [];
    keptlegs2 = 1:length(other.leg_names);
    keptlegs2(legpos2) = [];

    % We match the irreps according to the dependencies
    tmpirrepindices1 = 1:obj.irrep_number;
    tmpirrepindices2 = (1:other.irrep_number) + obj.irrep_number;
    matchlist_length = sum(cellfun(@(x) length(x), obj.dependencies(legpos1)));
    
    matchirrep_list1 = zeros(1,matchlist_length);
    matchirrep_list2 = zeros(1,matchlist_length);
    j = 1;
    for i = 1:length(legnames1)
        for depID = 1:length(obj.dependencies{legpos1(i)})
            irrepind1 = obj.dependencies{legpos1(i)}(depID);
            matchirrep_list1(j) = irrepind1;
            irrepind2 = other.dependencies{legpos2(i)}(depID);
            matchirrep_list2(j) = irrepind2;
            j = j + 1;
            if tmpirrepindices1(irrepind1) < tmpirrepindices2(irrepind2)
                % We change all the index values newirrepindices2(irrepind2) to newirrepindices1(irrpeind1)
                tmpirrepindices1(tmpirrepindices1 == tmpirrepindices2(irrepind2)) = tmpirrepindices1(irrepind1);
                tmpirrepindices2(tmpirrepindices2 == tmpirrepindices2(irrepind2)) = tmpirrepindices1(irrepind1);
            elseif tmpirrepindices1(irrepind1) > tmpirrepindices2(irrepind2)
                tmpirrepindices2(tmpirrepindices2 == tmpirrepindices1(irrepind1)) = tmpirrepindices2(irrepind2);
                tmpirrepindices1(tmpirrepindices1 == tmpirrepindices1(irrepind1)) = tmpirrepindices2(irrepind2);
            end
        end     
    end
    tmpirrepinds_all = unique([tmpirrepindices1,tmpirrepindices2]);

    % We sort the irrepindices to make the kept ones first and the
    % summed up after
    
    kept_tmpirrepinds = [];
    for i = 1:length(keptlegs1)
        kept_tmpirrepinds = [kept_tmpirrepinds,tmpirrepindices1(obj.dependencies{keptlegs1(i)})];
    end
    for i = 1:length(keptlegs2)
        kept_tmpirrepinds = [kept_tmpirrepinds,tmpirrepindices2(other.dependencies{keptlegs2(i)})];
    end
    kept_tmpirrepinds = unique(kept_tmpirrepinds,'sorted');
    sumup_tmpirrepinds = setdiff(tmpirrepinds_all,kept_tmpirrepinds,'sorted');
    neworder = [kept_tmpirrepinds,sumup_tmpirrepinds];
    new_irrep_number = length(kept_tmpirrepinds);   % the number of irreps in the product tensor
    
    %possibly slow below.
    new_irrepindices1 = tmpirrepindices1;
    new_irrepindices2 = tmpirrepindices2;
    for i = 1:length(neworder)
        new_irrepindices1(tmpirrepindices1 == neworder(i)) = i;
        new_irrepindices2(tmpirrepindices2 == neworder(i)) = i;
    end


    % Now we create the list of new leg names. If the name of a free (remaining) leg is listed in renaming{1} (legs
    % of obj), or renaming{2} (legs of other), then the new name is stored. Otherwise we store the name of the leg as
    % it is used in tensors obj or other. If there is ambiguity in the new leg names (the same name is defined for
    % two or more legs) error will be raised later.
    % At the same time we also create the list of new leg types.
    old_leg_names1 = obj.leg_names(keptlegs1);
    new_leg_names1 = old_leg_names1;
    for rnm = renaming{1}
        new_leg_names1{strcmp(rnm{1}{1},old_leg_names1)} = rnm{1}{2};
    end

    old_leg_names2 = other.leg_names(keptlegs2);
    new_leg_names2 = old_leg_names2;
    for rnm = renaming{2}
        new_leg_names2{strcmp(rnm{1}{1},old_leg_names2)} = rnm{1}{2};
    end

    new_leg_names = [new_leg_names1, new_leg_names2];
    new_leg_types = [obj.leg_types(keptlegs1), other.leg_types(keptlegs2)];

    if isempty(new_leg_names)
        just_a_number = true;
    else
        just_a_number = false;
    end

    % The dependencies of the new tensor are now determined. 
    if ~just_a_number
        new_dependencies1 = cell(1,length(keptlegs1));
        for i = 1:length(keptlegs1)
            new_dependencies1{i} = new_irrepindices1(obj.dependencies{keptlegs1(i)});
        end
        new_dependencies2 = cell(1,length(keptlegs2));
        for i = 1:length(keptlegs2)
            new_dependencies2{i} = new_irrepindices2(other.dependencies{keptlegs2(i)});
        end
        new_dependencies = [new_dependencies1, new_dependencies2];
    end


    % --- Output structure is built.
    out = struct();
    out.legpos1 = legpos1;
    out.legpos2 = legpos2;
    out.keptlegs1 = keptlegs1;
    out.keptlegs2 = keptlegs2;
    out.matchirrep_list1 = matchirrep_list1;
    out.matchirrep_list2 = matchirrep_list2;
    out.new_leg_names = new_leg_names;
    out.new_leg_types = new_leg_types;
    if ~just_a_number
        out.new_dependencies = new_dependencies;
    end
    out.new_irrepindices1 = new_irrepindices1;
    out.new_irrepindices2 = new_irrepindices2;
    out.new_irrep_number = new_irrep_number;
    out.just_a_number = just_a_number;

end

