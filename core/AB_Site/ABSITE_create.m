function obj = ABSITE_create(no_of_symmetries,sector_list,varargin)
            % Constructor of "Site" objects
            % symmetries:              cell array of Symmetry objects
            % sector_list:             cell array of sectors
            %                   format: {{[irrep labels1], taudim1},
            %                   {[irrep labels2], taudim2}, ...}
            
            parameters = struct();
            parameters.FERMIONS = false;
            parameters = parameter_updater(parameters,varargin);
            
            obj.type = 'AB_Site';
            obj.fermions = parameters.FERMIONS;
            obj.no_of_symmetries = no_of_symmetries;
            obj.sector_multiplets = containers.Map();
            for secID = 1:length(sector_list)
                if length(sector_list{secID}{1}) ~= obj.no_of_symmetries
                    error('Irrep labels must contain no_of_symmetries quantum numbers');
                end
                seckey = char(sector_list{secID}{1});
                obj.sector_multiplets(seckey) = sector_list{secID}{2};
            end
            obj.operators = containers.Map();
            obj.operator_irreps = containers.Map('UniformValues',false);
            obj.fermionicities = containers.Map();
            obj.operators('id') = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},obj.no_of_symmetries);
            for seckey_cell = obj.sector_multiplets.keys
                seckey = seckey_cell{1};
                obj.operators('id') = NTset_block(obj.operators('id'),{{'tau',1},{'tau~',1}},{double(seckey),double(seckey)},{'tau','tau~'},eye(obj.sector_multiplets(seckey)));
            end
            obj.fermionicities('id') = 1;
            obj.operator_irreps('id') = ones(1,no_of_symmetries);
            if obj.fermions        % We create a dummy phase operator (should be overwritten in general, but works for the auxiliary site of DMRG)
                obj.operators('ph') = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},obj.no_of_symmetries);
                for seckey_cell = obj.sector_multiplets.keys
                    seckey = seckey_cell{1};
                    obj.operators('ph') = NTset_block(obj.operators('ph'),{{'tau',1},{'tau~',1}},{double(seckey),double(seckey)},{'tau','tau~'},eye(obj.sector_multiplets(seckey)));
                end
            end
            obj.operator_irreps('ph') = ones(1,no_of_symmetries);
end

