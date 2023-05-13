function         [irrep_classes,...
                  irreps1_in_class,...
                  irreps2_in_class,...
                  just_a_number,...
                  legpos1,...
                  legpos2,...
                  match_classes,...
                  new_dependencies,...
                  new_leg_names,...
                  new_leg_types,...
                  sum_up_classes] = NTdot_init(obj,other,legnames1,legnames2,renaming)
% Performs preliminary book-keeping calculations for dot. Used also by LNAtensors.
            
            if length(legnames1) ~= length(legnames2)
                error('The number of multiplied legs are not the same in the two tensor!')
            end
            if length(unique(legnames1)) ~= length(legnames1) || ...
                   length(unique(legnames2)) ~= length(legnames2)
                error('Repetition in legnames1 or legnames2')
            end
            if obj.no_of_symmetries ~= other.no_of_symmetries
                error('The number of symmetries of the two tensors are not equal')
            end
            legpos1 = NTlocate_legs(obj, legnames1);
            legpos2 = NTlocate_legs(other, legnames2);
            for i = 1:length(legnames1)
                if obj.leg_types{legpos1(i)} == other.leg_types{legpos2(i)}
                    error('Legs with different types ("i","o") can only be contracted.')
                end
            end
            
            % First we perform the matching of irrep indices. We collect the irreps which are matched into lists
            % ("classes"). These classes corrsepond to the irrep indices of the resulting tensor. (Remark: some of them are
            % summed up later.) All classes are initialized with one member ({1,i} : i'th irrep of "self", {2,i} : i'th
            % irrep of "other"), and these classes are merged together, if there is a dependency match between them.
            irrep_classes = {};
            for i = 1:obj.irrep_number
                irrep_classes{end + 1} = {{1,i}};           %if it will be slow, we change it! (Copied from the python version.)
            end
            for i = 1:other.irrep_number
                irrep_classes{end + 1} = {{2,i}};
            end
            for i = 1:length(legnames1)
                if length(obj.dependencies{legpos1(i)}) ~= length(other.dependencies{legpos2(i)})
                    error('Dependency structures of the two tensors are not compatible.')
                end
                for j = 1:length(obj.dependencies{legpos1(i)})
                    element1 = {1, obj.dependencies{legpos1(i)}(j)};
                    element2 = {2, other.dependencies{legpos2(i)}(j)};       % WAS WRONG!
                    for classID = 1:length(irrep_classes)
                        if any(cellfun(@(x) isequal(x,element1),irrep_classes{classID}))
                            class1 = classID;
                        end
                        if any(cellfun(@(x) isequal(x,element2),irrep_classes{classID}))
                            class2 = classID;
                        end
                    end
                    if class1 > class2
                        irrep_classes{class2} = [irrep_classes{class2}, irrep_classes{class1}];
                        irrep_classes(class1) = [];
                    elseif class2 > class1
                        irrep_classes{class1} = [irrep_classes{class1}, irrep_classes{class2}];
                        irrep_classes(class2) = [];
                    end
                end
            end
            
            % Now we create the list of new leg names. If the name of a free (remaining) leg is listed in renaming{1} (legs
            % of obj), or renaming{2} (legs of other), then the new name is stored. Otherwise we store the name of the leg as
            % it is used in tensors obj or other. If there is ambiguity in the new leg names (the same name is defined for
            % two or more legs) error is raised.
            % At the same time we also create the list of new leg types.
            new_leg_names = cell(1,length(obj.leg_names) + length(other.leg_names) - 2*length(legnames1));
            new_leg_types = cell(1,length(obj.leg_names) + length(other.leg_names) - 2*length(legnames1));
            newind = 1;
            for ind = 1:length(obj.leg_names)
                if ~any(legpos1 == ind)
                    rnmfound = false;
                    name = obj.leg_names{ind};
                    for rnm = renaming{1}
                        if strcmp(rnm{1}{1},name) && ~rnmfound
                            name = rnm{1}{2};
                            rnmfound = true;
                        end
                    end
                    new_leg_types{newind} = obj.leg_types{ind};
                    new_leg_names{newind} = name;
                    newind = newind + 1;
                end
            end
            for ind = 1:length(other.leg_names)
                if ~any(legpos2 == ind)
                    rnmfound = false;
                    name = other.leg_names{ind};
                    for rnm = renaming{2}
                        if strcmp(rnm{1}{1},name) && ~rnmfound
                            name = rnm{1}{2};
                            rnmfound = true;
                        end
                    end
                    new_leg_types{newind} = other.leg_types{ind};
                    new_leg_names{newind} = name;
                    newind = newind + 1;
                end
            end
            if length(new_leg_names) ~= length(unique(new_leg_names))
                'new_leg_names = '
                new_leg_names
                error('Ambiguously defined leg names!')
            end
            if isempty(new_leg_names)
                just_a_number = true;
            else
                just_a_number = false;
            end
            
            % The dependencies of the new tensor are now determined. To do so, we search the dependency irreps within the
            % irrep_classes. The ID of the class will be the new irrep label in the new tensor.
            new_dependencies = cell(1,length(new_leg_names));
            newind = 1;
            for ind = 1:length(obj.leg_names)
                if ~any(legpos1 == ind)
                    deptmp = zeros(1, length(obj.dependencies{ind}));
                    for j = 1:length(obj.dependencies{ind})
                        elem = {1,obj.dependencies{ind}(j)};
                        for classID = 1:length(irrep_classes)
                            if any(cellfun(@(x) isequal(x,elem),irrep_classes{classID}))
                                deptmp(j) = classID;
                            end
                        end
                    end
                    new_dependencies{newind} = deptmp;
                    newind = newind + 1;
                end
            end
            for ind = 1:length(other.leg_names)
                if ~any(legpos2 == ind)
                    deptmp = zeros(1, length(other.dependencies{ind}));
                    for j = 1:length(other.dependencies{ind})
                        elem = {2, other.dependencies{ind}(j)};
                        for classID = 1:length(irrep_classes)
                            if any(cellfun(@(x) isequal(x,elem), irrep_classes{classID}))
                                deptmp(j) = classID;
                            end
                        end
                    end
                    new_dependencies{newind} = deptmp;
                    newind = newind + 1;
                end
            end
            
            % Now we search for irreps that should be summed up (irreps with no depending legs), and we also collect the
            % irrep classes where matching should performed, i.e. irrep classes with more than one member.
            sum_up_classes = [];
            match_classes = [];
            for classID = 1:length(irrep_classes)
                found = false;
                for deptmp = new_dependencies
                    if any(deptmp{1} == classID)
                        found = true;
                    end
                end
                if ~found
                    sum_up_classes(end+1) = classID;
                elseif length(irrep_classes{classID}) > 1
                    match_classes(end+1) = classID;
                end
            end
            for id = 1:length(new_dependencies)  % Decreasing the values of the remaining irrep labels to 1,2,3,..., i.e. irrep
                                                 % labels of the the sum-up classes are removed.
                for depID = 1:length(new_dependencies{id})
                    new_dependencies{id}(depID) = new_dependencies{id}(depID) ...
                        - sum(sum_up_classes < new_dependencies{id}(depID));
                end
            end

            % We slightly transform data in list irrep_classes. We collect irreps of obj (Irreps1_inClass) and irreps of
            % other (Irreps2_inClass) according to their classID.
            irreps1_in_class = cell(1,length(irrep_classes));      %MODIFIED 2021.06.17
            irreps2_in_class = cell(1,length(irrep_classes));      %MODIFIED 2021.06.17
            for classID = 1:length(irrep_classes)
                irreps1_in_class{classID} = [];
                irreps2_in_class{classID} = [];
                for elem = irrep_classes{classID}
                    if elem{1}{1} == 1
                        irreps1_in_class{classID}(end+1) = elem{1}{2};
                    else
                        irreps2_in_class{classID}(end+1) = elem{1}{2};
                    end
                end
                irreps1_in_class{classID} = sort(irreps1_in_class{classID});
                irreps2_in_class{classID} = sort(irreps2_in_class{classID});
            end

end

