function out = NTget_leg_sectors(obj,LegName)
    % Returns the sectors (list of quantum numbers, and the corresponding
    % dimensions) for a given leg
    % ---
    % obj:     the NAtensor structure
    % LegName  the name of the leg, whose sectors are listed
    % ---
    % out:     cell array (length: number of sectors) 
    %                - each cell contains a cell array of length 2:
    %                   * in the first cell the quantum numbers of the sector are given 
    %                        (cell array, length is the number of
    %                        dependencies of the leg)
    %                   * in the second cell the dimension of the sector is
    %                     given.
            legID = NTlocate_legs(obj,{LegName});
            key_layout = zeros(1,obj.no_of_symmetries*length(obj.dependencies{legID}));
            for depID = 1:length(obj.dependencies{legID})
                key_layout((obj.no_of_symmetries*(depID-1)+1):(obj.no_of_symmetries*depID)) = ...
                    (obj.no_of_symmetries*(obj.dependencies{legID}(depID)-1)+1):(obj.no_of_symmetries*obj.dependencies{legID}(depID));
            end
            outmap = containers.Map();
            for key = obj.data.keys
                seckey = key{1}(key_layout);
                if outmap.isKey(seckey)
                    tmpshape = size(obj.data(key{1}),legID);
                    if ~isequal(outmap(seckey), tmpshape)
                        error('NAtensor is corrupted. Different dimensions found at same leg sectors!')
                    end
                else
                    outmap(seckey) = size(obj.data(key{1}),legID);
                end
            end
            out = cell(1,outmap.Count);
            i = 1;
            for seckey = outmap.keys
                secirrep = cell(1,length(obj.dependencies{legID}));
                for depID = 1:length(obj.dependencies{legID})
                    secirrep{depID} = double(seckey{1}((obj.no_of_symmetries*(depID-1)+1):(obj.no_of_symmetries*depID)));
                end
                out{i}={secirrep,outmap(seckey{1})};
                i = i + 1;
            end
end

